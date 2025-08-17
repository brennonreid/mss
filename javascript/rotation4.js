/**
 * Rotate a 3D vector v by angle theta about axis k, using ONLY getUnitVectorFromAngle2 for cos/sin.
 * NOTE: k must be a UNIT vector. No per-frame normalization is performed here.
 * This is algebraically equivalent to the classic axis–angle (Rodrigues) apply-step,
 * but it does NOT use Math.sin/cos or any matrices/quaternions.
 *
 * @param {number[]} v - vector to rotate [vx, vy, vz]
 * @param {number[]} k - unit axis [kx, ky, kz]
 * @param {number} theta - angle in radians
 * @param {boolean} [precise=true] - pass-through to getUnitVectorFromAngle2 for Taylor depth
 * @returns {number[]} rotated vector [vx', vy', vz']
 */
function rotateAxisAngle(v, k, theta) {
  const { x: c, y: s } = getUnitVectorFromAngle2(theta, true); // c=cosθ, s=sinθ

  const kx = k[0], ky = k[1], kz = k[2];
  const vx = v[0], vy = v[1], vz = v[2];

  // k × v (cross)
  const cx = ky * vz - kz * vy;
  const cy = kz * vx - kx * vz;
  const cz = kx * vy - ky * vx;

  // k ⋅ v (dot)
  const kd = kx * vx + ky * vy + kz * vz;

  // v' = v*c + (k×v)*s + k*(k⋅v)*(1 - c)
  const oneMinusC = 1 - c;
  return [
    vx * c + cx * s + kx * kd * oneMinusC,
    vy * c + cy * s + ky * kd * oneMinusC,
    vz * c + cz * s + kz * kd * oneMinusC
  ];
}


function applyAxisRotateAxisAngle(x, y, z, axis, theta)
{
      const k = (axis === 'z') ? [0,0,1] : (axis === 'x') ? [1,0,0] : [0,1,0];
      const out = rotateAxisAngle([x, y, z], k, theta, true); // precise Taylor path in your lib
      return out; // [x', y', z']
}

function canonicalizeZXY(alpha, beta, gamma, precise = true) {
  const { x: cb } = getUnitVectorFromAngle2(beta, precise); // cos(beta)
  if (cb < 0) {
    alpha += Math.PI;
    beta   = Math.PI - beta;
    gamma += Math.PI;
  }
  return [alpha, beta, gamma];
}

// Apply yaw(Z), pitch(X), roll(Y) with your β-flip inside.
function applyYawPitchRoll(x, y, z, yaw, pitch, roll) {
  let alpha = yaw, beta = pitch, gamma = roll;
  [alpha, beta, gamma] = canonicalizeZXY(alpha, beta, gamma, true);
  let p = rotateAxisAngle([x, y, z], [0,0,1], alpha);
  p     = rotateAxisAngle(p,           [1,0,0], beta);
  p     = rotateAxisAngle(p,           [0,1,0], gamma);
  return p; // [x',y',z']
}

// Vector convenience
function applyYawPitchRollv(v, yaw, pitch, roll) {
  return applyYawPitchRoll(v[0], v[1], v[2], yaw, pitch, roll);
}

// Pose builder: pre-canonicalizes once; reuse apply() many times per frame.
function rotate(yaw, pitch, roll) {
  let alpha = yaw, beta = pitch, gamma = roll;
  [alpha, beta, gamma] = canonicalizeZXY(alpha, beta, gamma, true);
  return {
    apply(x, y, z) {
      let p = rotateAxisAngle([x, y, z], [0,0,1], alpha);
      p     = rotateAxisAngle(p,           [1,0,0], beta);
      p     = rotateAxisAngle(p,           [0,1,0], gamma);
      return p;
    },
    applyv(v) { return this.apply(v[0], v[1], v[2]); },
    angles() { return { yaw: alpha, pitch: beta, roll: gamma }; } // canonicalized
  };
}

// Precompute c/s once, then apply Z -> X -> Y as 2D rotates.
// Returns a pose object compatible with your draw code (has .apply and .angles()).
function makePoseZXY_CS(yaw, pitch, roll, precise = true) {
  let [alpha, beta, gamma] = canonicalizeZXY(yaw, pitch, roll, precise);

  const { x: cZ, y: sZ } = getUnitVectorFromAngle2(alpha, precise);
  const { x: cX, y: sX } = getUnitVectorFromAngle2(beta,  precise);
  const { x: cY, y: sY } = getUnitVectorFromAngle2(gamma, precise);

  return {
    apply(x, y, z) {
      // Z: rotate (x,y)
      const x1 = x * cZ - y * sZ;
      const y1 = x * sZ + y * cZ;
      const z1 = z;

      // X: rotate (y,z)
      const y2 = y1 * cX - z1 * sX;
      const z2 = y1 * sX + z1 * cX;
      const x2 = x1;

      // Y: rotate (x,z)
      const x3 = x2 * cY + z2 * sY;
      const z3 = -x2 * sY + z2 * cY;
      const y3 = y2;

      return [x3, y3, z3];
    },
    applyv(v) { return this.apply(v[0], v[1], v[2]); },
    angles() { return { yaw: alpha, pitch: beta, roll: gamma }; } // canonicalized
  };
}

// Precompute c/s once, then apply Z -> Y -> Z as 2D rotates.
// Uses your canonicalizeZXY treating (alpha=z1, beta=y, gamma=z2) — parity still correct.
function makePoseZYZ_CS(z1, y, z2, precise = true) {
  let [alpha, beta, gamma] = canonicalizeZXY(z1, y, z2, precise);

  const { x: cZ1, y: sZ1 } = getUnitVectorFromAngle2(alpha, precise);
  const { x: cY,  y: sY  } = getUnitVectorFromAngle2(beta,  precise);
  const { x: cZ2, y: sZ2 } = getUnitVectorFromAngle2(gamma, precise);

  return {
    apply(x, y, z) {
      // Z1: rotate (x,y)
      const x1 = x * cZ1 - y * sZ1;
      const y1 = x * sZ1 + y * cZ1;
      const z1 = z;

      // Y : rotate (x,z)
      const x2 = x1 * cY + z1 * sY;
      const z2 = -x1 * sY + z1 * cY;
      const y2 = y1;

      // Z2: rotate (x,y)
      const x3 = x2 * cZ2 - y2 * sZ2;
      const y3 = x2 * sZ2 + y2 * cZ2;
      const z3 = z2;

      return [x3, y3, z3];
    },
    applyv(v) { return this.apply(v[0], v[1], v[2]); },
    angles() { return { z1: alpha, y: beta, z2: gamma }; }
  };
}

