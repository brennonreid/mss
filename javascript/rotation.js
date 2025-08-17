// rotation_core.js — Matrix-free, trig-free rotation helpers (globals)
// Depends on customtrig.js providing:
//   TAU
//   getUnitVectorFromAngle2(theta[, precise])
//   atan2(y, x), asin(z)
// cos/sin via your oracle (tiny gated renorm to stay unit)
function cosSinFromAngle(theta){
  const u = getUnitVectorFromAngle2(theta, /*precise=*/false); // {x,y}
  let c = u.x, s = u.y;
  const m2 = c*c + s*s;
  if (Math.abs(m2 - 1.0) > 1e-12) {
    const inv = 1.0 / Math.sqrt(m2);
    c *= inv; s *= inv;
  }
  return { c, s };
}

// ---- Direct axis–angle rotate (Rodrigues) ----
// axis can be non-unit; we normalize once inside.
// ---- Direct axis–angle rotate (Rodrigues) ----
// axis can be non-unit; normalize only if needed.
function rotateAxisAngle(kx, ky, kz, theta, x, y, z){
  const kn2 = kx*kx + ky*ky + kz*kz;             // faster than hypot + avoids sqrt when unit
  if (Math.abs(kn2 - 1.0) > 1e-12) {
    const inv = 1 / Math.sqrt(kn2);
    kx *= inv; ky *= inv; kz *= inv;
  }

  const uv = getUnitVectorFromAngle2(theta, /*precise=*/false); // {x: cosθ, y: sinθ}
  const c = uv.x, s = uv.y;

  // k × v
  const cx = ky * z - kz * y;
  const cy = kz * x - kx * z;
  const cz = kx * y - ky * x;

  // k · v
  const kd = kx * x + ky * y + kz * z;

  const oneMinusC = 1 - c;
  return [
    x * c + cx * s + kx * kd * oneMinusC,
    y * c + cy * s + ky * kd * oneMinusC,
    z * c + cz * s + kz * kd * oneMinusC,
  ];
}


// Axis helpers (unit axes by definition)
function rotateX(theta, x, y, z){ return rotateAxisAngle(1,0,0, theta, x,y,z); }
function rotateY(theta, x, y, z){ return rotateAxisAngle(0,1,0, theta, x,y,z); }
function rotateZ(theta, x, y, z){ return rotateAxisAngle(0,0,1, theta, x,y,z); }

// Chained convenience (your demo order Z→X→Y)
function rotateZXY(alpha, beta, gamma, x, y, z){
  let v = rotateZ(alpha, x,y,z);
  v = rotateX(beta, v[0], v[1], v[2]);
  v = rotateY(gamma, v[0], v[1], v[2]);
  return v;
}

// ---- Euler readout (ZYX) directly from axis–angle (no matrix object) ----
function eulerZYXFromAxisAngle(kx, ky, kz, theta){
  // normalize axis
  const kn = Math.hypot(kx, ky, kz) || 1;
  kx /= kn; ky /= kn; kz /= kn;

  const { c, s } = cosSinFromAngle(theta);
  const C = 1 - c;

  // Rodrigues components we need (not building/storing a matrix)
  const r00 = c + kx*kx*C;
  const r01 = kx*ky*C - kz*s;
  const r02 = kx*kz*C + ky*s;

  const r10 = ky*kx*C + kz*s;
  const r11 = c + ky*ky*C;
  const r12 = ky*kz*C - kx*s;

  const r20 = kz*kx*C - ky*s;
  const r21 = kz*ky*C + kx*s;
  const r22 = c + kz*kz*C;

  const pitch = asin(-r20);
  const cp = Math.hypot(r00, r10);

  let yaw, roll;
  if (cp < 1e-6){
    yaw  = atan2(-r01, r11);
    roll = 0;
  } else {
    yaw  = atan2(r10, r00);
    roll = atan2(r21, r22);
  }
  return { yaw, pitch, roll };
}

// ---- Compose two axis–angle rotations (no matrices/quats exposed) ----
// Returns equivalent axis–angle of Rb ∘ Ra.
// Uses half-angle representation internally with sqrt + atan2 (allowed).
function composeAxisAngle(ax, ay, az, aTheta, bx, by, bz, bTheta){
  // normalize axes
  let an = Math.hypot(ax, ay, az) || 1; ax/=an; ay/=an; az/=an;
  let bn = Math.hypot(bx, by, bz) || 1; bx/=bn; by/=bn; bz/=bn;

  // bring angles to (-π, π] for stable half-angle sign
  const PI = TAU * 0.5;
  function wrapPi(t){ t = (t + PI) % TAU; return (t <= 0) ? t + PI : t - PI; }
  aTheta = wrapPi(aTheta);
  bTheta = wrapPi(bTheta);

  // cos/sin(θ) via oracle → half-angles by identities:
  // cos(θ/2)=sqrt((1+cosθ)/2), sin(θ/2)=sign(θ)*sqrt((1-cosθ)/2)
  const ca = getUnitVectorFromAngle2(aTheta, false).x;
  const cb = getUnitVectorFromAngle2(bTheta, false).x;

  const cah = Math.sqrt(Math.max(0, 0.5*(1 + ca)));
  let  sah = Math.sqrt(Math.max(0, 0.5*(1 - ca))); if (aTheta < 0) sah = -sah;

  const cbh = Math.sqrt(Math.max(0, 0.5*(1 + cb)));
  let  sbh = Math.sqrt(Math.max(0, 0.5*(1 - cb))); if (bTheta < 0) sbh = -sbh;

  // "quaternion-like" multiply (w, v) without exposing quats
  // qa = (cah, sah * â), qb = (cbh, sbh * b̂)  → qc = qb * qa
  const vax = ax * sah, vay = ay * sah, vaz = az * sah;
  const vbx = bx * sbh, vby = by * sbh, vbz = bz * sbh;

  const w  = cbh*cah - (vbx*vax + vby*vay + vbz*vaz);
  const vx = cbh*vax + cah*vbx + (vay*vbz - vaz*vby);
  const vy = cbh*vay + cah*vby + (vaz*vbx - vax*vbz);
  const vz = cbh*vaz + cah*vbz + (vax*vby - vay*vbx);

  const sn = Math.hypot(vx, vy, vz);
  if (sn < 1e-18) {
    // nearly identity or 180° about opposite axes; pick a stable axis
    return { kx: 1, ky: 0, kz: 0, theta: 0 };
  }

  // θ = 2 * atan2(|v|, w), k̂ = v / |v|
  const theta = 2 * atan2(sn, w);
  const inv = 1 / sn;
  return { kx: vx*inv, ky: vy*inv, kz: vz*inv, theta };
}

// ---- rotate by Z, X, Y angles without matrices (convenience) ----
function rotateZXY(alpha, beta, gamma, x, y, z){
  let v = rotateZ(alpha, x,y,z);
  v = rotateX(beta, v[0], v[1], v[2]);
  v = rotateY(gamma, v[0], v[1], v[2]);
  return v;
}

// (optional) rotation-vector helpers — still matrix-free
function rotateVectorFromRotVec(rx, ry, rz, x, y, z){
  const theta = Math.hypot(rx, ry, rz);
  if (theta === 0) return [x, y, z];
  const inv = 1 / theta;
  return rotateAxisAngle(rx*inv, ry*inv, rz*inv, theta, x, y, z);
}

function eulerZYXFromRotVec(rx, ry, rz){
  const theta = Math.hypot(rx, ry, rz);
  if (theta === 0) return { yaw:0, pitch:0, roll:0 };
  const inv = 1 / theta;
  return eulerZYXFromAxisAngle(rx*inv, ry*inv, rz*inv, theta);
}