# Modular Symbolic System 

This is a trigonometry-free, drift-free, rotation system built for high-performance applications.

Faster than trig, matrix, and quaternion systems (see demos).

---

## What This Is

This project offers a new approach to rotation that eliminates common pain points found in traditional systems.

By combining:

- **Trig-free computation** via Symbolically Enriched Lookup Tables (SE-LUT)
- **3D rotation** via phase-tracked axis–angle logic
- **Deterministic control** via per-axis symbolic stepping

…it delivers a **drift-free**, **gimbal-free**, **high-performance** rotation engine, faster and more predictable than trig, matrix, or quaternion systems.

It contains a minimal set of JS/Java/C++ demos for:

- **Trig-computation** via Symbolically Enriched Lookup Tables (SE-LUT)
- **2D rotation** via complex multiply
- **3D rotation** via axis-angle (Rodrigues) logic
- **Per-axis phase control** (no global coupling)

---

## What Makes This Trig System Special?

This project does **not** use a traditional trigonometric lookup table.

Instead, it uses a **Symbolically Enriched Lookup Table (SE-LUT)**, where each anchor:

- Encodes full phase and vector context (`θ`, `cosθ`, `sinθ`)
- Supports symbolic Taylor expansion from that anchor
- Enables phase-tracked, drift-free trigonometric evaluation

As a result, even with just **128 anchors**, the trig system achieves:

- **Machine-precision output** (RMS ≈ 2.6e-16)
- **30%+ speedup** over native `Math.sin/cos`

---

## Rotation System (Rodrigues + Phase)

On top of the SE-LUT trigonometry layer, this project builds a full 3D rotation engine using:

- **Rodrigues-style axis–angle rotation**
- **Per-axis phase tracking** (no accumulation, no coupling)
- **Pure unit vector operations** (no quats, no matrices)

**Why Rodrigues?**

Rodrigues rotation offers precise 3D vector rotation using only vector algebra, making it a perfect match for phase-tracked systems. By separating axis and angle as symbolic phase values, this method avoids the pitfalls of quaternions (like normalization) while remaining compact and intuitive.

Key properties:

- **Drift-free**: integer phase stepping, no loss over time  
- **Matrix-free**: no persistent 3×3 state  
- **Quaternion-free**: no need for 4D logic or normalization  
- **Gimbal-free**: no Euler collapse or ambiguity  
- **Fast**: beats quaternions and matrices in apply speed  

---

## Benchmarks: Trig vs Native, Rotation vs Quaternion

<details>
<summary><strong>Click to view full results</strong></summary>

### Trig System: SE-LUT vs Native `Math.cos/sin`

```text
--- UNIFORM (10M samples) ---
Native Math.cos/sin        : 918.60 ms (10.89 M/s)
Custom (MSS precise=true)  : 642.70 ms (15.56 M/s)

--- SMALL-DELTA (STEP=2π/512) ---
Native Math.cos/sin        : 923.90 ms
Custom (MSS precise=true)  : 668.30 ms

RMS Error: 2.577e−16
Max Error: 9.105e−16
Speedup: ~1.4× faster

=== 2,000,000 runs ===
Z-only     : MSS = 112.7 ms  | Quat = 186.7 ms  → 1.66× faster
Y-only     : MSS = 105.7 ms  | Quat = 190.1 ms  → 1.80× faster
ZXY combo  : MSS = 461.7 ms  | Quat = 567.7 ms  → 1.23× faster

Max component deviation: ≤ 2.50e−15
Checksums match perfectly.
```
</details>

---

## License & Patent Notice (Interim)

This project is under **patent pending** status.

This repository is provided for **non-commercial evaluation and research only**. A formal license will follow.

### Permitted Uses

- Clone/view/run demos  
- Non-commercial use only  
- Attribution required: `brennonreid/mss`  

### Prohibited Uses

- Commercial products or services  
- Redistribution or sublicensing  
- Derivative works  
- Use in ML training or datasets  

---

**AS IS** – No warranty.  
© 2025 Brennon Reid. All rights reserved.
