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
- **Quaternion-free:** no quaternion math required for composition.
- **Gimbal-free:** no Euler-angle degeneracy in operation.
- **Fast:** beats native trig and quaternion apply in our benches.

## Accuracy & Speed
We’re faster and just as accurate.

- **Unit vectors (cos/sin):** ~**1.4×** faster on average. Accuracy: RMS **2.577e−16**, max **9.105e−16** vs `Math.cos/sin`.
- **3D rotations (vs quaternions):** up to **1.8×** faster. Accuracy: checksums match; max component deviation ≤ **2.50e−15**.
<details>
<summary><strong>Full benchmark data</strong></summary>

#### Unit vectors
```text
--- UNIFORM (10,000,000 samples) ---
Native Math.cos/sin        : 918.60 ms   (10.89 M/s)   checksum=2977.770125814061
Custom (MSS precise=true)  : 642.70 ms   (15.56 M/s)   checksum=2977.770125814480

--- SMALL-DELTA (STEP=2π/512) (10,000,000 samples) ---
Native Math.cos/sin        : 923.90 ms   (10.82 M/s)   checksum=937.344599554581
Custom (MSS precise=true)  : 668.30 ms   (14.96 M/s)   checksum=937.344599554971

RMS / maxAbs on 2,000,000 uniform samples → rms=2.577e-16  max=9.105e-16
Throughput gain: UNIFORM 1.43× (−30.0% time); SMALL-DELTA 1.38× (−27.7% time).

3D rotations
text
Copy
Edit
=== Grouped run  N=2,000,000  seed=1337 ===
Z-only         MSS(dir)=112.700 ms  ns/op=56.350  Quat=186.700 ms  ns/op=93.350  chkD=3.321e+5  chkQ=3.321e+5  maxDev=2.498e-15
Y-only         MSS(dir)=105.700 ms  ns/op=52.850  Quat=190.100 ms  ns/op=95.050  chkD=-6.659e+5 chkQ=-6.659e+5 maxDev=1.443e-15
ZXY composite  MSS(dir)=461.700 ms  ns/op=230.850 Quat=567.700 ms  ns/op=283.850 chkD=-430.608 chkQ=-430.608 maxDev=2.442e-15
--- Summary ---
Z-only          MSS(dir)=112.700 ms  Quat=186.700 ms  Δ=74.000 ms
Y-only          MSS(dir)=105.700 ms  Quat=190.100 ms  Δ=84.400 ms
ZXY             MSS(dir)=461.700 ms  Quat=567.700 ms  Δ=106.000 ms
Throughput gain: Z-only 1.66× (−39.6% time), Y-only 1.80× (−44.4%), ZXY 1.23× (−18.7%).
```
</details> 

## License & Patent Notice (Interim)

**Status: Patent Pending**  
This repository is currently for **non-commercial evaluation and research only**. No use is permitted in commercial products or services. A permanent license will be added later.

### Permitted Uses
- View, clone, compile locally, and run the included demos.
- Use is restricted to **non-commercial evaluation and research only**.
- Attribution is required: `brennonreid/mss`.

### Prohibited Uses
- **Commercial use:** No incorporation into products or services.
- **Redistribution or hosting:** No public hosting, sublicensing, or redistribution.
- **Derivatives:** No creation or distribution of derivative works.
- **ML/AI:** No inclusion in datasets or ML training.

**Formal Notice:**  
No express or implied rights are granted under any patents or patent applications related to Modular Symbolic System (MSS)/RSR or related work.

Provided “AS IS”, without warranties or conditions of any kind. Use at your own risk. This interim notice controls until a formal `LICENSE` file is added.

© 2025 Brennon Reid. All rights reserved.