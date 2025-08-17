const TAU = 2 * Math.PI;
const QUADTAU = TAU / 4; 
const NUM_ANCHORS_QUADRANT = 128;
const DEG_STEP_QUADRANT = QUADTAU / NUM_ANCHORS_QUADRANT;
const NUM_ANCHORS_FULL = 4 * NUM_ANCHORS_QUADRANT;
const DEG_STEP_FULL = TAU / NUM_ANCHORS_FULL;
const STEP = TAU / NUM_ANCHORS_FULL;

const quadrantTransforms = [
    (x, y) => ({ x, y }),
    (x, y) => ({ x: -y, y: x }),
    (x, y) => ({ x: -x, y: -y }),
    (x, y) => ({ x: y, y: -x }),
];

// Apply the quadrant transformation based on the given angle
function applyQuadrantTransform(x, y, quadrant) {
    return quadrantTransforms[quadrant](x, y);
}

// Taylor Series Expansion for Cosine
function taylorCosFloat(x, depth = 2) {
    let term = 1, sum = 1;
    for (let i = 1; i <= depth; i++) {
        term *= -x * x / ((2 * i - 1) * (2 * i));
        sum += term;
    }
    return sum;
}

// Taylor Series Expansion for Sine
function taylorSinFloat(x, depth = 2) {
    let term = x, sum = x;
    for (let i = 1; i <= depth; i++) {
        term *= -x * x / ((2 * i) * (2 * i + 1));
        sum += term;
    }
    return sum;
}

// Create Canonical Anchors (for base angle steps)
const canonicalAnchors = [];
for (let i = 0; i < NUM_ANCHORS_QUADRANT; i++) {
    const angle = i * DEG_STEP_QUADRANT;
    canonicalAnchors.push({
        theta: angle,
        cos: taylorCosFloat(angle, 55),
        sin: taylorSinFloat(angle, 55),
    });
}

// Build Full Circle Anchors from Canonical Anchors
const fullCircleAnchors = [];
for (let i = 0; i < NUM_ANCHORS_FULL; i++) {
    const angle = i * DEG_STEP_FULL;
    const quadrant = Math.floor(angle / QUADTAU);
    const theta = angle - quadrant * QUADTAU;
    const canonicalIndex = Math.round(theta / DEG_STEP_QUADRANT);
    const canonicalAnchor = canonicalAnchors[Math.min(canonicalIndex, NUM_ANCHORS_QUADRANT - 1)];
    const { x, y } = applyQuadrantTransform(canonicalAnchor.cos, canonicalAnchor.sin, quadrant);
    fullCircleAnchors.push({ 
        theta: angle, 
        cos: x, 
        sin: y 
    });
}

// Custom Sine Function Using Taylor Expansion
function getCustomSin(angle, precise = false) {
    angle %= TAU;
    if (angle < 0) angle += TAU;

    const anchorIndex = Math.floor(angle / STEP);
    const anchor = fullCircleAnchors[anchorIndex];
    const delta = angle - anchorIndex * STEP;
    const δ = delta, δ2 = δ * δ, δ3 = δ * δ2, δ5 = δ3 * δ2, δ7 = δ5 * δ2;

    const sinValue = precise
        ? δ - δ3 * (1 / 6 - δ2 * (1 / 120 - δ2 * (1 / 5040 - δ2 / 362880)))
        : δ - δ3 * (1 / 6 - δ2 * (1 / 120 - δ2 / 5040));

    return sinValue;
}

// Custom Cosine Function Using Taylor Expansion
function getCustomCos(angle, precise = false) {
    angle %= TAU;
    if (angle < 0) angle += TAU;

    const anchorIndex = Math.floor(angle / STEP);
    const anchor = fullCircleAnchors[anchorIndex];
    const delta = angle - anchorIndex * STEP;
    const δ = delta, δ2 = δ * δ, δ4 = δ2 * δ2, δ6 = δ4 * δ2;

    const cosValue = precise
        ? 1 - δ2 * (1 / 2 - δ2 * (1 / 24 - δ2 * (1 / 720 - δ2 / 40320)))
        : 1 - δ2 * (1 / 2 - δ2 * (1 / 24 - δ2 / 720));

    return cosValue;
}

// Unit Vector Calculation with Precise Taylor Approximation
function getUnitVectorFromAngle2(angle, precise = false) {
    angle %= TAU;
    if (angle < 0) angle += TAU;

    const anchorIndex = Math.floor(angle / STEP);
    const anchor = fullCircleAnchors[anchorIndex];
    const delta = angle - anchorIndex * STEP;
    const δ = delta, δ2 = δ * δ;

	// inline taylor for speed
    const dx = precise
        ? 1 - δ2 * (0.5 - δ2 * (1 / 24 - δ2 * (1 / 720 - δ2 / 40320)))
        : 1 - δ2 * (0.5 - δ2 * (1 / 24 - δ2 / 720));

    const dy = precise
        ? δ * (1 - δ2 * (1 / 6 - δ2 * (1 / 120 - δ2 / 5040)))
        : δ * (1 - δ2 * (1 / 6 - δ2 / 120));

    const x = dx * anchor.cos - dy * anchor.sin;
    const y = dx * anchor.sin + dy * anchor.cos;

	return { x, y };
}

  function wrapTau(a){ a %= TAU; return a < 0 ? a + TAU : a; }

  // Trig-free angle from (x,y). Returns angle in (-π, π].
  // M: coarse samples, refine: bisection steps. Tweak if desired.
  function angle_from_xy_symbolic(x, y, M = 64, refine = 18){
    if (x === 0 && y === 0) return 0;

    // normalize v
    const n = Math.hypot(x, y) || 1;
    const vx = x / n, vy = y / n;

    // coarse: choose best sample by dot
    let bestIdx = 0, bestDot = -Infinity;
    for (let i = 0; i < M; i++){
      const a = (i * TAU) / M;
      const u = getUnitVectorFromAngle2(a, /*precise=*/false);
      const inv = 1 / (Math.hypot(u.x, u.y) || 1);
      const dot = (u.x * inv) * vx + (u.y * inv) * vy;
      if (dot > bestDot){ bestDot = dot; bestIdx = i; }
    }

    // bracket: neighbors of best sample
    let lo = ((bestIdx - 1 + M) % M) * TAU / M;
    let hi = ((bestIdx + 1) % M) * TAU / M;
    if (hi < lo) hi += TAU;

    // refine by bisection on cross sign (u × v)
    for (let k = 0; k < refine; k++){
      const mid = 0.5 * (lo + hi);
      const a = wrapTau(mid);
      const u = getUnitVectorFromAngle2(a, /*precise=*/false);
      const inv = 1 / (Math.hypot(u.x, u.y) || 1);
      const ux = u.x * inv, uy = u.y * inv;
      const cross = ux * vy - uy * vx; // >0 ⇒ u is CCW of v
      if (cross > 0) lo = mid; else hi = mid;
    }

    let ang = wrapTau(0.5 * (lo + hi));
    if (ang > TAU * 0.5) ang -= TAU;  // map to (-π, π]
    return ang;
  }

  // Trig-free atan2/asin built on the symbolic angle finder.
  function atan2(y, x){ return angle_from_xy_symbolic(x, y); }

  function asin(z){
    const HALFPI = TAU * 0.25;
    if (z >=  1) return  HALFPI;
    if (z <= -1) return -HALFPI;
    const c = Math.sqrt(Math.max(0, 1 - z*z)); // |c| = sqrt(1 - z^2)
    return atan2(z, c);
  }

