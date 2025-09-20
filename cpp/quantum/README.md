# btrig / quantum

Tiny, header-only helpers that turn **fixed-point phases** (Q64.64 / Q128.128 “turns”) into **quantum gate programs**.  
Focus: **exact dyadic angles**, clean **OpenQASM-style emit**, and fast **batch phase tooling** that pairs with [`btrig`](../btrig).

> Why: many quantum workflows are “phase-heavy.” If your classical side already stores angles as fractions of a turn, you can emit gates **exactly** (as `2*pi*k/2^n`) or **numerically** with a single switch—no ambiguous rounding.

---

## Features

- **Fixed-point phases**  
  - `Phase64` (Q64.64) and `Phase128` (Q128.128) with add/half/neg and degree/turn/radian converters.
- **Angle modes (one API, two encodings)**  
  - `AngleEmit::Dyadic` → prints `rz(2*pi*k/2^n) q[i];` (no float rounding).  
  - `AngleEmit::Float64` → prints numeric radians (double).
- **OpenQASM-style emitter**  
  - Minimal `OpenQASM3Emitter`: `qreg`, `rz/ry/rx`, `barrier`, `measure`, string/file dump.
  - Raw expression helpers: `rz_expr/rx_expr/ry_expr("2*pi*(...)", q[i])`.
- **Schedules & helpers**  
  - `make_phase_qft(nbits)` → descending dyadic phases for QFT-like patterns.  
  - `make_phase_linear(base, step, count)` → affine phase sweeps.  
- **Optional fast DDS streamer** (define `BQ_ENABLE_DDS`)  
  - `dds_stream_q64(phases, callback)` batches through your SIMD kernel.

This folder is deliberately **small**: it’s not a simulator or transpiler. It’s a **precise, readable phase-to-gates bridge**.

---

## Quick start

```cpp
#include "btrig/btrig.h"          // Q64.64 math (turns)
#include "quantum/btrig_quantum.h"

using namespace bq; // namespace of the quantum helpers

int main() {
    OpenQASM3Emitter qasm;

    auto q = qasm.qreg("q", 2);

    // 45 degrees as a Q64.64 phase (one turn = 360°)
    Phase64 p = Phase64::from_degrees_ld(45.0L);

    // 1) Exact dyadic: prints rz(2*pi*k/2^64) q[0];
    emit_rz(qasm, q[0], p, AngleEmit::Dyadic);

    // 2) Numeric radians (double)
    emit_rx(qasm, q[1], p, AngleEmit::Float64);

    // Barrier + measurement for completeness
    qasm.barrier(q);
    qasm.measure(q, "c");

    // Get the program as a string (or qasm.dump_to_file("prog.qasm");)
    std::string text = qasm.str();
    // ... use with your favorite toolchain / simulator
}
Build (header-only)

bash
Copy code
g++ -O3 -std=c++17 -I ../../cpp demo.cpp -o demo
# (Optional) If you enable DDS streaming that calls SIMD:
# g++ -O3 -std=c++17 -mavx2 -mfma -DBQ_ENABLE_DDS -I ../../cpp demo.cpp -o demo
API snapshot
Phases
using Phase64 = PhaseQ<64>;

using Phase128 = PhaseQ<128>;

Constructors:

Phase64::from_degrees_ld(long double deg);

Phase64::from_turns_q64(btrig::uq64 turns);

Phase64::from_radians_ld(long double rad);

Ops:

p.add(q), p.neg(), p.half() (adds ½ turn), p.to_radians_ld(), p.to_turns_ld()

Emitters
OpenQASM3Emitter:

QReg qreg(std::string name, int n);

void barrier(const QReg&);

void measure(const QReg&, std::string cbitBase);

std::string str() const;

void dump_to_file(const char* path) const;

Gate emission (normalized phases):

emit_rz(OpenQASM3Emitter&, QRef, Phase64, AngleEmit mode);

emit_rx(...), emit_ry(...)

Exact helpers: emit_rz_dyadic(...), emit_rx_dyadic(...), emit_ry_dyadic(...)

Raw expression (when you want full control):

rz_expr(OpenQASM3Emitter&, QRef, const std::string& expr);

Same for rx_expr/ry_expr.

Schedules
std::vector<Phase64> make_phase_qft(int nbits);

std::vector<Phase64> make_phase_linear(Phase64 base, Phase64 step, int count);

Optional: DDS streaming (batch phase → sin/cos)
Define BQ_ENABLE_DDS and include your SIMD kernel (btrig_simd.h).

dds_stream_q64(phases, [&](double c, double s, int i){ ... });

Exact vs numeric angles
Dyadic (“exact”): angles are emitted as 2*pi*k/2^n. If your control stack or simulator treats this symbolically, you avoid any float rounding (nice for QFT, phase kickback, and provable identities).

Numeric: prints a double radians literal. Good for general use and tools that require concrete numbers.

Pick per-gate with AngleEmit::Dyadic or AngleEmit::Float64.

Precision notes
Internally, phases are turns: 1.0 turn = τ radians.

Phase64 (Q64.64) gives sub-attosecond-like granularity in turns (≈ 5.4e-20).

When emitting Dyadic, the printed form is exact (symbolic).

When emitting Float64, the literal is double—use Dyadic if you need exactness.

Minimal CMake (header-only)
Add this in your project’s CMake and include the repo path:

cmake
Copy code
# btrig core (turns + math)
add_library(btrig INTERFACE)
target_include_directories(btrig INTERFACE ${CMAKE_CURRENT_LIST_DIR}/..)
target_compile_features(btrig INTERFACE cxx_std_17)

# quantum helpers
add_library(btrig_quantum INTERFACE)
target_link_libraries(btrig_quantum INTERFACE btrig)
target_include_directories(btrig_quantum INTERFACE ${CMAKE_CURRENT_LIST_DIR}/..)

# optional SIMD if you enable DDS streaming
add_library(btrig_simd INTERFACE)
target_link_libraries(btrig_simd INTERFACE btrig)
target_include_directories(btrig_simd INTERFACE ${CMAKE_CURRENT_LIST_DIR}/..)
target_compile_options(btrig_simd INTERFACE -mavx2 -mfma)
Roadmap (non-blocking)
Reversible integer add / add-constant emitters (Cuccaro-style).

More schedule builders (chirps, windowed ramps).

QASM printer toggles (strict QASM 2 vs QASM 3 flavor).

Error-budget helpers for mapping fixed-point phase to hardware dials.

Contributing
Issues and PRs welcome. Please include a short example (input phases → expected gate text) in new tests.

See DEDICATION.md for the homage behind this work.
