// btrig.js

// --------------------
// Configuration
// --------------------
const TRIG_ANCHORS_BASE_POW = 5;
const NUM_ANCHORS_QUADRANT = 1 << TRIG_ANCHORS_BASE_POW;
const NUM_ANCHORS_FULL = NUM_ANCHORS_QUADRANT * 4;
const ANCHOR_MASK = NUM_ANCHORS_FULL - 1;

// --------------------
// Core constants
// --------------------
const TAU = 6.283185307179586476925286766559; // 2*pi
const QUADTAU = TAU / 4.0;
const STEP = TAU / NUM_ANCHORS_FULL;
const ONE_OVER_STEP = 1.0 / STEP;

const DEG_STEP_QUADRANT = QUADTAU / NUM_ANCHORS_QUADRANT;
const DEG_STEP_FULL = TAU / NUM_ANCHORS_FULL;

// Precomputed inverse factorials for interpolation
const INV_2 = 1.0 / 2.0;
const INV_6 = 1.0 / 6.0;
const INV_24 = 1.0 / 24.0;
const INV_120 = 1.0 / 120.0;
const INV_720 = 1.0 / 720.0;
const INV_5040 = 1.0 / 5040.0;
const INV_40320 = 1.0 / 40320.0;
const INV_362880 = 2.755731922398589065255731922e-06;  // 1/9!
const INV_3628800 = 2.755731922398589065255731922e-07; // 1/10!

const INV_TAU = 1.0 / TAU;

// Anchor array: flat [theta, cos, sin] triplets
// length = NUM_ANCHORS_FULL * 3
const fullCircleAnchors = new Float64Array(NUM_ANCHORS_FULL * 3);

// --------------------
// Internal helpers
// --------------------
function applyQuadrantTransform(inx, iny, quadrant) {
  switch (quadrant & 3) {
    case 0: return { x: inx,  y: iny  };
    case 1: return { x: -iny, y: inx  };
    case 2: return { x: -inx, y: -iny };
    default:return { x: iny,  y: -inx };
  }
}

function wrapTau(a) {
  // C++ version computes floor-like q carefully. JS version is simple and not hot-path.
  const k = a * INV_TAU;
  let q = Math.trunc(k);
  if (a < 0.0 && q !== k) q -= 1;
  return a - q * TAU;
}

function taylorCosFloat(x, depth) {
  let term = 1.0;
  let sum = 1.0;
  for (let i = 1; i <= depth; ++i) {
    term *= -x * x / ((2.0 * i - 1.0) * (2.0 * i));
    sum += term;
  }
  return sum;
}

function taylorSinFloat(x, depth) {
  let term = x;
  let sum = x;
  for (let i = 1; i <= depth; ++i) {
    term *= -x * x / ((2.0 * i) * (2.0 * i + 1.0));
    sum += term;
  }
  return sum;
}

function taylorSinCosFloat2(angle, depth) {
  const a = wrapTau(angle);
  const quadrant = Math.floor(a / QUADTAU);
  const x = a - quadrant * QUADTAU;

  const sinVal = taylorSinFloat(x, depth);
  const cosVal = taylorCosFloat(x, depth);

  const t = applyQuadrantTransform(cosVal, sinVal, quadrant);
  return { sin: t.y, cos: t.x };
}

function initFullCircleAnchors() {
  for (let i = 0; i < NUM_ANCHORS_FULL; ++i) {
    const angle = i * DEG_STEP_FULL;

    const quadrant = Math.floor(angle / QUADTAU);
    const x = angle - quadrant * QUADTAU;

    const cosVal = taylorCosFloat(x, 55);
    const sinVal = taylorSinFloat(x, 55);

    const t = applyQuadrantTransform(cosVal, sinVal, quadrant);

    const base = i * 3;
    fullCircleAnchors[base + 0] = angle;
    fullCircleAnchors[base + 1] = t.x; // cos
    fullCircleAnchors[base + 2] = t.y; // sin
  }
}

// Eager init to match static C++ initializer
initFullCircleAnchors();

// Careful floor toward -inf using trunc fix, matching C++ pattern:
// int i = (int)t; i -= (t < 0.0);
function floorTowardNegInfFromTruncLike(t) {
  let i = Math.trunc(t);
  if (t < 0.0 && i !== t) i -= 1;
  return i;
}

// --------------------
// Public API (mirrors your C++ names)
// --------------------
const btrig = {
  // constants (exported for parity/debug)
  TRIG_ANCHORS_BASE_POW,
  NUM_ANCHORS_QUADRANT,
  NUM_ANCHORS_FULL,
  ANCHOR_MASK,
  TAU,
  QUADTAU,
  STEP,
  ONE_OVER_STEP,
  DEG_STEP_QUADRANT,
  DEG_STEP_FULL,
  INV_2,
  INV_6,
  INV_24,
  INV_120,
  INV_720,
  INV_5040,
  INV_40320,
  INV_362880,
  INV_3628800,
  INV_TAU,

  // expose anchors if you want to inspect
  fullCircleAnchors,

  applyQuadrantTransform,
  wrapTau,
  taylorCosFloat,
  taylorSinFloat,
  taylorSinCosFloat2,

  sin: function(angle, precise, outObj) {
    // outObj: { value: number } to emulate pass-by-ref
    const t = angle * ONE_OVER_STEP;

    let i = floorTowardNegInfFromTruncLike(t);
    // mask assumes power-of-two table size; use modulo if you prefer:
    const idx = i & ANCHOR_MASK;

    const d = (t - i) * STEP;
    const d2 = d * d;

    const base = idx * 3;
    const a1 = fullCircleAnchors[base + 1]; // cos
    const a2 = fullCircleAnchors[base + 2]; // sin

    let dx, dy;
    if (precise) {
      dx = 1.0 - d2 * (INV_2 -
           d2 * (INV_24 -
           d2 * (INV_720 -
           d2 *  INV_40320)));
      dy = d * (1.0 -
           d2 * (INV_6 -
           d2 * (INV_120 -
           d2 *  INV_5040)));
    } else {
      dx = 1.0 - d2 * INV_2;
      dy = d;
    }
    outObj.value = dx * a2 + dy * a1;
  },

  cos: function(angle, precise, outObj) {
    const t = angle * ONE_OVER_STEP;

    let i = floorTowardNegInfFromTruncLike(t);
    const idx = i & ANCHOR_MASK;

    const d = (t - i) * STEP;
    const d2 = d * d;

    const base = idx * 3;
    const a1 = fullCircleAnchors[base + 1]; // cos
    const a2 = fullCircleAnchors[base + 2]; // sin

    let dx, dy;
    if (precise) {
      dx = 1.0 - d2 * (INV_2 -
           d2 * (INV_24 -
           d2 * (INV_720 -
           d2 *  INV_40320)));
      dy = d * (1.0 -
           d2 * (INV_6 -
           d2 * (INV_120 -
           d2 *  INV_5040)));
    } else {
      dx = 1.0 - d2 * INV_2;
      dy = d;
    }
    outObj.value = dx * a1 - dy * a2;
  },

  sincos: function(angle, precise, outCosObj, outSinObj) {
    const t = angle * ONE_OVER_STEP;

    let i = floorTowardNegInfFromTruncLike(t);
    const idx = i & ANCHOR_MASK;

    const d = (t - i) * STEP;
    const d2 = d * d;

    const base = idx * 3;
    const a1 = fullCircleAnchors[base + 1]; // cos
    const a2 = fullCircleAnchors[base + 2]; // sin

    let dx, dy;
    if (precise) {
      dx = 1.0 - d2 * (INV_2 -
           d2 * (INV_24 -
           d2 * (INV_720 -
           d2 *  INV_40320)));
      dy = d * (1.0 -
           d2 * (INV_6 -
           d2 * (INV_120 -
           d2 *  INV_5040)));
    } else {
      dx = 1.0 - d2 * INV_2;
      dy = d;
    }

    outCosObj.value = dx * a1 - dy * a2;
    outSinObj.value = dx * a2 + dy * a1;
  },

  atan2: function(y, x) {
    if (x === 0.0 && y === 0.0) return 0.0;

    const ax = (y < 0.0 ? -y : y);
    const ay = (x < 0.0 ? -x : x);

    const swap = (ax > ay);
    const num = swap ? ay : ax;
    const den = swap ? ax : ay;
    const z = num / den;

    const z2 = z * z;
    const a = z * (0.9998660 +
            z2 * (-0.3302995 +
            z2 * (0.1801410 +
            z2 * (-0.0851330 +
                   0.0208351 * z2)))));

    const angle = swap ? (1.5707963267948966 - a) : a;

    if (x >= 0.0) {
      return (y >= 0.0) ? angle : -angle;
    } else {
      return (y >= 0.0) ? (3.141592653589793 - angle) : (angle - 3.141592653589793);
    }
  },

  // Assumes IEEE-754 doubles like the C++ comment.
  fast_sqrt: function(x) {
    if (x <= 0.0) return 0.0;
    // JS has no bit-cast trick; emulate fast inverse sqrt via Newton steps on 1/sqrt
    // Start with 1/Math.sqrt(x) as a decent guess; if you truly want no Math.sqrt,
    // you can seed with a rough power approximation.
    let y = 1.0 / Math.sqrt(x);
    const xhalf = 0.5 * x;
    // Two Newton refinements on inverse sqrt
    y = y * (1.5 - xhalf * y * y);
    y = y * (1.5 - xhalf * y * y);
    return x * y;
  },

  asin: function(z) {
    const HALFPI = TAU * 0.25;
    if (z >= 1.0)  return HALFPI;
    if (z <= -1.0) return -HALFPI;

    const t = 1.0 - z * z;
    const c = (t > 0.0) ? btrig.fast_sqrt(t) : 0.0;
    return btrig.atan2(z, c);
  },

  acos: function(x) {
    if (x >= 1.0)  return 0.0;
    if (x <= -1.0) return TAU * 0.5;

    const t = 1.0 - x * x;
    const s = (t > 0.0) ? btrig.fast_sqrt(t) : 0.0;
    return btrig.atan2(s, x);
  },

  tan: function(angle) {
    const c = { value: 0.0 };
    const s = { value: 0.0 };
    btrig.sincos(angle, true, c, s);
    return s.value / c.value;
  },

  atan: function(y) {
    return btrig.atan2(y, 1.0);
  }
};

// CommonJS and ESM friendly exports
// Node (CommonJS):
if (typeof module !== "undefined" && module.exports) {
  module.exports = btrig;
}
// ESM default export:
export default btrig;
