Deterministic, phase-driven rotation.
**Drift-free. Matrix-free. Quaternion-free. Gimbal-free.**
Faster in trig and in rotations (see demos/benches).

## What this is
A minimal set of JS/Java/C++ demos and helpers for reliable 2D/3D rotation using:
- Phase + LUT unit vectors (no runtime `Math.sin/Math.cos` in the hot path)
- 2D complex multiply for planar rotation
- Axis-angle (Rodrigues) for 3D rotation
- Per-axis phase (no global reset coupling)

## Why it matters
- **Drift-free:** deterministic phase stepping; no accumulation bleed.
- **Matrix-free:** no persistent 3×3 stacks.
- **Quaternion-free:** no q math required for composition.
- **Gimbal-free:** no Euler-angle degeneracy in operation.
- **Fast:** beats native trig and quaternion apply in our benches.
