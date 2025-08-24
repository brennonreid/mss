// rotation.js — Matrix-free helpers wired to btrig2.js
// Depends on a global `btrig` object providing:
//   btrig.TAU
//   btrig.sincos(theta, precise) -> {cos, sin}
//   btrig.atan2(y, x)
//   btrig.asin(z)

// Shorthands to mirror the old names
const TAU = btrig.TAU;
const atan2 = (y, x) => btrig.atan2(y, x);
const asin  = (z)    => btrig.asin(z);

// --- Compact inline cos/sin oracle (unit-safe) ---
function cosSinFromAngle(theta, precise = false) {
  const r = btrig.sincos(theta, precise);
  return { c: r.cos, s: r.sin };
}

// --- Direct axis–angle rotation (Rodrigues) ---
function rotateAxisAngle(kx, ky, kz, theta, x, y, z) {
  const r = btrig.sincos(theta, /*precise=*/false);
  const c = r.cos, s = r.sin;

  const cx = ky * z - kz * y;
  const cy = kz * x - kx * z;
  const cz = kx * y - ky * x;

  const kd = kx * x + ky * y + kz * z;
  const oneMinusC = 1 - c;

  return [
    x * c + cx * s + kx * kd * oneMinusC,
    y * c + cy * s + ky * kd * oneMinusC,
    z * c + cz * s + kz * kd * oneMinusC,
  ];
}

// --- Axis helpers (unit axes only) ---
function rotateX(theta, x, y, z) { return rotateAxisAngle(1, 0, 0, theta, x, y, z); }
function rotateY(theta, x, y, z) { return rotateAxisAngle(0, 1, 0, theta, x, y, z); }
function rotateZ(theta, x, y, z) { return rotateAxisAngle(0, 0, 1, theta, x, y, z); }

// --- Chained rotate (Z → X → Y) ---
function rotateZXY(alpha, beta, gamma, x, y, z) {
  let v = rotateZ(alpha, x, y, z);
  v = rotateX(beta, v[0], v[1], v[2]);
  v = rotateY(gamma, v[0], v[1], v[2]);
  return v;
}

// --- Euler readout (ZYX) from axis–angle ---
function eulerZYXFromAxisAngle(kx, ky, kz, theta) {
  const { c, s } = cosSinFromAngle(theta);
  const C = 1 - c;

  const r00 = c + kx * kx * C;
  const r01 = kx * ky * C - kz * s;
  const r02 = kx * kz * C + ky * s;

  const r10 = ky * kx * C + kz * s;
  const r11 = c + ky * ky * C;
  const r12 = ky * kz * C - kx * s;

  const r20 = kz * kx * C - ky * s;
  const r21 = kz * ky * C + kx * s;
  const r22 = c + kz * kz * C;

  const pitch = asin(-r20);
  const cp = Math.hypot(r00, r10);

  let yaw, roll;
  if (cp < 1e-6) {
    yaw  = atan2(-r01, r11);
    roll = 0;
  } else {
    yaw  = atan2(r10, r00);
    roll = atan2(r21, r22);
  }

  return { yaw, pitch, roll };
}

// --- Compose axis–angle → axis–angle (no matrix/quat exposed) ---
function composeAxisAngle(ax, ay, az, aTheta, bx, by, bz, bTheta) {
  function wrapPi(t) {
    t = (t + TAU / 2) % TAU;
    return (t <= 0) ? t + TAU / 2 : t - TAU / 2;
  }

  aTheta = wrapPi(aTheta);
  bTheta = wrapPi(bTheta);

  // cos(a), cos(b) via btrig
  const ca = btrig.sincos(aTheta, false).cos;
  const cb = btrig.sincos(bTheta, false).cos;

  // half-angle cos/sin from full-angle cos (sign from original angle)
  const cah = Math.sqrt(Math.max(0, 0.5 * (1 + ca)));
  let sah = Math.sqrt(Math.max(0, 0.5 * (1 - ca)));
  if (aTheta < 0) sah = -sah;

  const cbh = Math.sqrt(Math.max(0, 0.5 * (1 + cb)));
  let sbh = Math.sqrt(Math.max(0, 0.5 * (1 - cb)));
  if (bTheta < 0) sbh = -sbh;

  // Quaternion-like composition in minimal form
  const vax = ax * sah, vay = ay * sah, vaz = az * sah;
  const vbx = bx * sbh, vby = by * sbh, vbz = bz * sbh;

  const w  = cbh * cah - (vbx * vax + vby * vay + vbz * vaz);
  const vx = cbh * vax + cah * vbx + (vay * vbz - vaz * vby);
  const vy = cbh * vay + cah * vby + (vaz * vbx - vax * vbz);
  const vz = cbh * vaz + cah * vbz + (vax * vby - vay * vbx);

  const sn = Math.hypot(vx, vy, vz);
  if (sn < 1e-18) return { kx: 1, ky: 0, kz: 0, theta: 0 };

  const theta = 2 * atan2(sn, w);
  const inv = 1 / sn;
  return { kx: vx * inv, ky: vy * inv, kz: vz * inv, theta };
}

// --- Rotate via rotation vector ---
function rotateVectorFromRotVec(rx, ry, rz, x, y, z) {
  const theta = Math.hypot(rx, ry, rz);
  if (theta === 0) return [x, y, z];
  const inv = 1 / theta;
  return rotateAxisAngle(rx * inv, ry * inv, rz * inv, theta, x, y, z);
}

function eulerZYXFromRotVec(rx, ry, rz) {
  const theta = Math.hypot(rx, ry, rz);
  if (theta === 0) return { yaw: 0, pitch: 0, roll: 0 };
  const inv = 1 / theta;
  return eulerZYXFromAxisAngle(rx * inv, ry * inv, rz * inv, theta);
}
