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

## Temporary Evaluation License & Patent Notice (Interim)

**Status:** Patent pending. A permanent license will be added later. Until then, these terms apply to the entire repository.

- **Permitted:** View, clone, compile locally, and run the included demos **for non-commercial evaluation and research only**, with attribution to `brennonreid/mss`.
- **Not permitted:** Redistribution, public hosting, commercial use, sublicensing, incorporation into products/services, inclusion in datasets or ML training, or creation/distribution of derivative works **without prior written permission**.
- **No patent license:** No express or implied rights are granted under any patents or patent applications related to Modular Symbolic System (MSS)/RSR or related work.
- **No contributions:** Pull requests are not being accepted at this time. Bug reports via Issues are welcome.
- **No warranty:** Provided “AS IS”, without warranties or conditions of any kind. Use at your own risk.
- **Supersedes:** This interim notice controls until a formal `LICENSE` is added.

© 2025 Brennon Reid. All rights reserved.
