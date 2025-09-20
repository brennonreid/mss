#pragma once
#ifndef BTRIG_SIMD_H
#define BTRIG_SIMD_H

// SIMD extensions for btrig.h using AVX2 (+optional FMA) on Linux/GCC/Clang.
// - Vectorizes the residual polynomial + anchor rotation for batches of phases.
// - Keeps the original scalar Q64.64 phase splitter (idx + rem) per element
//   for correctness; computes the heavy math (dx,dy and final rotate) in SIMD.
// - Uses double precision for SIMD lanes for speed; anchors are converted once
//   from btrig::fullCircleAnchors (long double) to double.
// - Falls back to scalar loops if AVX2 is unavailable.
//
// Build example (Linux):
//   g++ -O3 -march=native -mavx2 -mfma demo.cpp -o demo
// or explicitly:
//   g++ -O3 -mavx2 -mfma demo.cpp -o demo
//
// Usage:
//   #include "btrig.h"
//   #include "btrig_simd.h"
//   using namespace btrig;
//   std::vector<uq64> phases(n);
//   std::vector<double> c(n), s(n);
//   simd::sincos_q64(phases.data(), n, /*precise=*/true, c.data(), s.data());
//
// Notes:
// * AVX2 path processes 4 doubles per iteration. Tail uses a 2-lane SSE2 block
//   (with AVX2 gathers) and then scalar.
// * Precision: SIMD uses double-precision (53-bit) anchors + polynomials for
//   speed. The scalar path in btrig.h remains long double. If you require
//   long-double outputs for specific lanes, call the scalar API for those.
// * Alignment: anchors are 32-byte aligned. Output pointers need not be
//   aligned; we choose aligned stores when the pointer is 32-byte aligned.
//   For peak throughput, pass 32-byte aligned cosOut/sinOut.
// * FMA: compile-time gated via __FMA__. Runtime has_fma() is provided only
//   for information; executing FMA instructions requires compiling with -mfma
//   (or -march=native where available).

#include <immintrin.h>
#include <cstdint>
#include <cstring>
#include <cstddef>  // alignment checks
#include <vector>
#include <cmath>
#include <mutex>

#include "btrig.h"

namespace btrig { namespace simd {

// -----------------------------------------------------------------------------
// Internal: double-precision copies of anchors for AVX2 gathers
// -----------------------------------------------------------------------------
alignas(32) static double _anchorsCosD[NUM_ANCHORS_FULL];
alignas(32) static double _anchorsSinD[NUM_ANCHORS_FULL];

static inline void _ensure_anchorsD_init() {
    static bool inited = false;
    if (inited) return;
    for (int i = 0; i < NUM_ANCHORS_FULL; ++i) {
        _anchorsCosD[i] = static_cast<double>(fullCircleAnchors[i][1]);
        _anchorsSinD[i] = static_cast<double>(fullCircleAnchors[i][2]);
    }
    inited = true;
}

// -----------------------------------------------------------------------------
// Runtime feature check (GCC/Clang):
// -----------------------------------------------------------------------------
static inline bool has_avx2() {
#if defined(__GNUC__)
    return __builtin_cpu_supports("avx2");
#else
    return false;
#endif
}
static inline bool has_fma() {
#if defined(__GNUC__)
    return __builtin_cpu_supports("fma");
#else
    return false;
#endif
}

// Pointer alignment helper
static inline bool is_aligned_32(const void* p) {
    return (reinterpret_cast<uintptr_t>(p) & 31u) == 0u;
}

// -----------------------------------------------------------------------------
// SIMD kernels (AVX2): process 4 lanes of double
// -----------------------------------------------------------------------------
namespace detail {

// Convert 4 signed 64-bit residuals to double radians via K_RESIDUAL
static inline __m256d load_d_from_srem4(const int64_t srem[4]) {
    const double K = static_cast<double>(btrig::K_RESIDUAL);
    return _mm256_set_pd((double)srem[3] * K,
                         (double)srem[2] * K,
                         (double)srem[1] * K,
                         (double)srem[0] * K);
}

// Gather anchors (cos,sin) using 32-bit indices
static inline void gather_anchor4(const int idx2[4], __m256d& cA, __m256d& sA) {
    __m128i vi = _mm_loadu_si128(reinterpret_cast<const __m128i*>(idx2)); // 4x int32
    cA = _mm256_i32gather_pd(_anchorsCosD, vi, 8); // 8-byte stride
    sA = _mm256_i32gather_pd(_anchorsSinD, vi, 8);
}

// Compute residual micro-rotation (dx,dy) for vector d.
// precise=true : cos(d) to 8!, sin(d) to 7!
// precise=false: cos(d) 4th,  sin(d) 3rd
static inline void residual4(__m256d d, bool precise, __m256d& dx, __m256d& dy) {
    __m256d z = _mm256_mul_pd(d, d);

    if (precise) {
        // tc = -1/8!; tc = tc*z + 1/6!; tc = tc*z - 1/4!; tc = tc*z + 1/2!
        __m256d tc = _mm256_set1_pd(-1.0/40320.0);
#if defined(__FMA__)
        tc = _mm256_fmadd_pd(tc, z, _mm256_set1_pd(1.0/720.0));
        tc = _mm256_fmadd_pd(tc, z, _mm256_set1_pd(-1.0/24.0));
        tc = _mm256_fmadd_pd(tc, z, _mm256_set1_pd(1.0/2.0));
        dx = _mm256_fnmadd_pd(z, tc, _mm256_set1_pd(1.0)); // 1 - z*tc
#else
        tc = _mm256_add_pd(_mm256_mul_pd(tc, z), _mm256_set1_pd(1.0/720.0));
        tc = _mm256_add_pd(_mm256_mul_pd(tc, z), _mm256_set1_pd(-1.0/24.0));
        tc = _mm256_add_pd(_mm256_mul_pd(tc, z), _mm256_set1_pd(1.0/2.0));
        dx = _mm256_sub_pd(_mm256_set1_pd(1.0), _mm256_mul_pd(z, tc));
#endif

        // ts = -1/7!; ts = ts*z + 1/5!; ts = ts*z - 1/3!; ts = ts*z + 1
        __m256d ts = _mm256_set1_pd(-1.0/5040.0);
#if defined(__FMA__)
        ts = _mm256_fmadd_pd(ts, z, _mm256_set1_pd(1.0/120.0));
        ts = _mm256_fmadd_pd(ts, z, _mm256_set1_pd(-1.0/6.0));
        ts = _mm256_fmadd_pd(ts, z, _mm256_set1_pd(1.0));
        dy = _mm256_mul_pd(d, ts);
#else
        ts = _mm256_add_pd(_mm256_mul_pd(ts, z), _mm256_set1_pd(1.0/120.0));
        ts = _mm256_add_pd(_mm256_mul_pd(ts, z), _mm256_set1_pd(-1.0/6.0));
        ts = _mm256_add_pd(ts, _mm256_set1_pd(1.0));
        dy = _mm256_mul_pd(d, ts);
#endif
    } else {
        // fast: cos(d) ~ 1 - d^2/2 + d^4/24 ; sin(d) ~ d - d^3/6
        __m256d tc = _mm256_set1_pd(-1.0/24.0);
#if defined(__FMA__)
        tc = _mm256_fmadd_pd(tc, z, _mm256_set1_pd(1.0/2.0));
        dx = _mm256_fnmadd_pd(z, tc, _mm256_set1_pd(1.0));
#else
        tc = _mm256_add_pd(_mm256_mul_pd(tc, z), _mm256_set1_pd(1.0/2.0));
        dx = _mm256_sub_pd(_mm256_set1_pd(1.0), _mm256_mul_pd(z, tc));
#endif

        __m256d ts = _mm256_set1_pd(-1.0/6.0);
#if defined(__FMA__)
        ts = _mm256_fmadd_pd(ts, z, _mm256_set1_pd(1.0));
        dy = _mm256_mul_pd(d, ts);
#else
        ts = _mm256_add_pd(_mm256_mul_pd(ts, z), _mm256_set1_pd(1.0));
        dy = _mm256_mul_pd(d, ts);
#endif
    }
}

// Final rotate: (c,s) = (dx*cA - dy*sA, dx*sA + dy*cA)
static inline void rotate4(__m256d dx, __m256d dy, __m256d cA, __m256d sA,
                           __m256d& c, __m256d& s) {
#if defined(__FMA__)
    __m256d neg_dy = _mm256_sub_pd(_mm256_set1_pd(0.0), dy);
    c = _mm256_fmadd_pd(dx, cA, _mm256_mul_pd(neg_dy, sA));
    s = _mm256_fmadd_pd(dx, sA, _mm256_mul_pd(dy, cA));
#else
    c = _mm256_sub_pd(_mm256_mul_pd(dx, cA), _mm256_mul_pd(dy, sA));
    s = _mm256_add_pd(_mm256_mul_pd(dx, sA), _mm256_mul_pd(dy, cA));
#endif
}

} // namespace detail

// -----------------------------------------------------------------------------
// Public AVX2 batch APIs
// -----------------------------------------------------------------------------
static inline void sincos_q64_avx2(const uq64* phase_q64,
                                   int n,
                                   bool precise,
                                   double* cosOut,
                                   double* sinOut) {
    _ensure_anchorsD_init();

#if defined(__AVX2__)
    int i = 0;
    for (; i + 3 < n; i += 4) {
        // Scalar index+remainder + right-edge pivot per lane
        int idx2[4];
        int64_t srem[4];
        for (int k = 0; k < 4; ++k) {
            int idx; uint64_t rem;
            btrig::phaseQ64_to_idx_rem(phase_q64[i + k], idx, rem);
            uint64_t sign = rem >> 63; // top bit
            idx2[k] = (idx + (int)sign) & btrig::ANCHOR_MASK;
            srem[k] = (int64_t)rem; // signed remainder around pivot
        }

        // d (radians), gather anchors, residual, rotate
        __m256d d  = detail::load_d_from_srem4(srem);
        __m256d cA, sA; detail::gather_anchor4(idx2, cA, sA);
        __m256d dx, dy; detail::residual4(d, precise, dx, dy);
        __m256d c, s;   detail::rotate4(dx, dy, cA, sA, c, s);

        // store
        if (is_aligned_32(&cosOut[i])) _mm256_store_pd(&cosOut[i], c); else _mm256_storeu_pd(&cosOut[i], c);
        if (is_aligned_32(&sinOut[i])) _mm256_store_pd(&sinOut[i], s); else _mm256_storeu_pd(&sinOut[i], s);
    }

    // partial tail: 2-lane SSE2 with AVX2 gathers, then scalar
    for (; i + 1 < n; i += 2) {
        int idx2[2];
        int64_t srem[2];
        for (int k = 0; k < 2; ++k) {
            int idx; uint64_t rem;
            btrig::phaseQ64_to_idx_rem(phase_q64[i + k], idx, rem);
            uint64_t sign = rem >> 63;
            idx2[k] = (idx + (int)sign) & btrig::ANCHOR_MASK;
            srem[k] = (int64_t)rem;
        }
        const double K = static_cast<double>(btrig::K_RESIDUAL);
        __m128d d = _mm_set_pd((double)srem[1] * K, (double)srem[0] * K);
        __m128d z = _mm_mul_pd(d, d);
        __m128i vi = _mm_loadu_si128(reinterpret_cast<const __m128i*>(idx2));
        __m128d cA = _mm_i32gather_pd(_anchorsCosD, vi, 8);
        __m128d sA = _mm_i32gather_pd(_anchorsSinD, vi, 8);
        __m128d dx, dy;
        if (precise) {
#if defined(__FMA__)
            __m128d tc = _mm_set1_pd(-1.0/40320.0);
            tc = _mm_fmadd_pd(tc, z, _mm_set1_pd(1.0/720.0));
            tc = _mm_fmadd_pd(tc, z, _mm_set1_pd(-1.0/24.0));
            tc = _mm_fmadd_pd(tc, z, _mm_set1_pd(1.0/2.0));
            dx = _mm_fnmadd_pd(z, tc, _mm_set1_pd(1.0));
            __m128d ts = _mm_set1_pd(-1.0/5040.0);
            ts = _mm_fmadd_pd(ts, z, _mm_set1_pd(1.0/120.0));
            ts = _mm_fmadd_pd(ts, z, _mm_set1_pd(-1.0/6.0));
            ts = _mm_fmadd_pd(ts, z, _mm_set1_pd(1.0));
            dy = _mm_mul_pd(d, ts);
#else
            __m128d tc = _mm_set1_pd(-1.0/40320.0);
            tc = _mm_add_pd(_mm_mul_pd(tc, z), _mm_set1_pd(1.0/720.0));
            tc = _mm_add_pd(_mm_mul_pd(tc, z), _mm_set1_pd(-1.0/24.0));
            tc = _mm_add_pd(tc, _mm_set1_pd(1.0/2.0));
            dx = _mm_sub_pd(_mm_set1_pd(1.0), _mm_mul_pd(z, tc));
            __m128d ts = _mm_set1_pd(-1.0/5040.0);
            ts = _mm_add_pd(_mm_mul_pd(ts, z), _mm_set1_pd(1.0/120.0));
            ts = _mm_add_pd(_mm_mul_pd(ts, z), _mm_set1_pd(-1.0/6.0));
            ts = _mm_add_pd(ts, _mm_set1_pd(1.0));
            dy = _mm_mul_pd(d, ts);
#endif
        } else {
#if defined(__FMA__)
            __m128d tc = _mm_set1_pd(-1.0/24.0);
            tc = _mm_fmadd_pd(tc, z, _mm_set1_pd(1.0/2.0));
            dx = _mm_fnmadd_pd(z, tc, _mm_set1_pd(1.0));
            __m128d ts = _mm_set1_pd(-1.0/6.0);
            ts = _mm_fmadd_pd(ts, z, _mm_set1_pd(1.0));
            dy = _mm_mul_pd(d, ts);
#else
            __m128d tc = _mm_set1_pd(-1.0/24.0);
            tc = _mm_add_pd(_mm_mul_pd(tc, z), _mm_set1_pd(1.0/2.0));
            dx = _mm_sub_pd(_mm_set1_pd(1.0), _mm_mul_pd(z, tc));
            __m128d ts = _mm_set1_pd(-1.0/6.0);
            ts = _mm_add_pd(_mm_mul_pd(ts, z), _mm_set1_pd(1.0));
            dy = _mm_mul_pd(d, ts);
#endif
        }
        __m128d c = _mm_sub_pd(_mm_mul_pd(dx, cA), _mm_mul_pd(dy, sA));
        __m128d s = _mm_add_pd(_mm_mul_pd(dx, sA), _mm_mul_pd(dy, cA));
        _mm_storeu_pd(&cosOut[i], c);
        _mm_storeu_pd(&sinOut[i], s);
    }
    for (; i < n; ++i) {
        long double c, s;
        btrig::sincos_q64(phase_q64[i], precise, c, s);
        cosOut[i] = (double)c;
        sinOut[i] = (double)s;
    }
#else
    // No AVX2 at compile-time: fallback to scalar implementation
    for (int i = 0; i < n; ++i) {
        long double c, s;
        btrig::sincos_q64(phase_q64[i], precise, c, s);
        cosOut[i] = (double)c;
        sinOut[i] = (double)s;
    }
#endif
}

// Convenience wrappers: cosine-only / sine-only
static inline void cos_q64_avx2(const uq64* phase_q64, int n, bool precise, double* cosOut) {
#if defined(__AVX2__)
    _ensure_anchorsD_init();
    int i = 0;
    for (; i + 3 < n; i += 4) {
        int idx2[4];
        int64_t srem[4];
        for (int k = 0; k < 4; ++k) {
            int idx; uint64_t rem;
            btrig::phaseQ64_to_idx_rem(phase_q64[i + k], idx, rem);
            uint64_t sign = rem >> 63;
            idx2[k] = (idx + (int)sign) & btrig::ANCHOR_MASK;
            srem[k] = (int64_t)rem;
        }
        __m256d d  = detail::load_d_from_srem4(srem);
        __m256d cA, sA; detail::gather_anchor4(idx2, cA, sA);
        __m256d dx, dy; detail::residual4(d, precise, dx, dy);
        __m256d c, s;  detail::rotate4(dx, dy, cA, sA, c, s);
        if (is_aligned_32(&cosOut[i])) _mm256_store_pd(&cosOut[i], c); else _mm256_storeu_pd(&cosOut[i], c);
    }
    // 2-lane tail with gathers
    for (; i + 1 < n; i += 2) {
        int idx2[2];
        int64_t srem[2];
        for (int k = 0; k < 2; ++k) {
            int idx; uint64_t rem;
            btrig::phaseQ64_to_idx_rem(phase_q64[i + k], idx, rem);
            uint64_t sign = rem >> 63;
            idx2[k] = (idx + (int)sign) & btrig::ANCHOR_MASK;
            srem[k] = (int64_t)rem;
        }
        const double K = (double)btrig::K_RESIDUAL;
        __m128d d = _mm_set_pd((double)srem[1]*K, (double)srem[0]*K);
        __m128d z = _mm_mul_pd(d, d);
        __m128i vi = _mm_loadu_si128(reinterpret_cast<const __m128i*>(idx2));
        __m128d cA = _mm_i32gather_pd(_anchorsCosD, vi, 8);
        __m128d sA = _mm_i32gather_pd(_anchorsSinD, vi, 8);
        __m128d dx, dy;
        if (precise) {
#if defined(__FMA__)
            __m128d tc=_mm_set1_pd(-1.0/40320.0);
            tc=_mm_fmadd_pd(tc, z, _mm_set1_pd(1.0/720.0));
            tc=_mm_fmadd_pd(tc, z, _mm_set1_pd(-1.0/24.0));
            tc=_mm_fmadd_pd(tc, z, _mm_set1_pd(1.0/2.0));
            dx=_mm_fnmadd_pd(z, tc, _mm_set1_pd(1.0));
            __m128d ts=_mm_set1_pd(-1.0/5040.0);
            ts=_mm_fmadd_pd(ts, z, _mm_set1_pd(1.0/120.0));
            ts=_mm_fmadd_pd(ts, z, _mm_set1_pd(-1.0/6.0));
            ts=_mm_fmadd_pd(ts, z, _mm_set1_pd(1.0));
            dy=_mm_mul_pd(d, ts);
#else
            __m128d tc=_mm_set1_pd(-1.0/40320.0);
            tc=_mm_add_pd(_mm_mul_pd(tc, z), _mm_set1_pd(1.0/720.0));
            tc=_mm_add_pd(_mm_mul_pd(tc, z), _mm_set1_pd(-1.0/24.0));
            tc=_mm_add_pd(tc, _mm_set1_pd(1.0/2.0));
            dx=_mm_sub_pd(_mm_set1_pd(1.0), _mm_mul_pd(z, tc));
            __m128d ts=_mm_set1_pd(-1.0/5040.0);
            ts=_mm_add_pd(_mm_mul_pd(ts, z), _mm_set1_pd(1.0/120.0));
            ts=_mm_add_pd(_mm_mul_pd(ts, z), _mm_set1_pd(-1.0/6.0));
            ts=_mm_add_pd(ts, _mm_set1_pd(1.0));
            dy=_mm_mul_pd(d, ts);
#endif
        } else {
#if defined(__FMA__)
            __m128d tc=_mm_set1_pd(-1.0/24.0);
            tc=_mm_fmadd_pd(tc, z, _mm_set1_pd(1.0/2.0));
            dx=_mm_fnmadd_pd(z, tc, _mm_set1_pd(1.0));
            __m128d ts=_mm_set1_pd(-1.0/6.0);
            ts=_mm_fmadd_pd(ts, z, _mm_set1_pd(1.0));
            dy=_mm_mul_pd(d, ts);
#else
            __m128d tc=_mm_set1_pd(-1.0/24.0);
            tc=_mm_add_pd(_mm_mul_pd(tc, z), _mm_set1_pd(1.0/2.0));
            dx=_mm_sub_pd(_mm_set1_pd(1.0), _mm_mul_pd(z, tc));
            __m128d ts=_mm_set1_pd(-1.0/6.0);
            ts=_mm_add_pd(_mm_mul_pd(ts, z), _mm_set1_pd(1.0));
            dy=_mm_mul_pd(d, ts);
#endif
        }
        __m128d c=_mm_sub_pd(_mm_mul_pd(dx, cA), _mm_mul_pd(dy, sA));
        _mm_storeu_pd(&cosOut[i], c);
    }
    for (; i < n; ++i) {
        long double c, s; btrig::sincos_q64(phase_q64[i], precise, c, s); cosOut[i]=(double)c; }
#else
    for (int i=0;i<n;++i){ long double c,s; btrig::sincos_q64(phase_q64[i], precise, c, s); cosOut[i]=(double)c; }
#endif
}

static inline void sin_q64_avx2(const uq64* phase_q64, int n, bool precise, double* sinOut) {
#if defined(__AVX2__)
    _ensure_anchorsD_init();

    int i = 0;

    // 4-lane AVX2 main loop
    for (; i + 3 < n; i += 4) {
        int     idx2[4];
        int64_t srem[4];

        // scalar phase split + nearest-anchor pivot per lane
        for (int k = 0; k < 4; ++k) {
            int idx; uint64_t rem;
            btrig::phaseQ64_to_idx_rem(phase_q64[i + k], idx, rem);
            const uint64_t sign = rem >> 63;
            idx2[k]  = (idx + (int)sign) & btrig::ANCHOR_MASK;
            srem[k]  = (int64_t)rem;
        }

        __m256d d  = detail::load_d_from_srem4(srem);
        __m256d cA, sA; detail::gather_anchor4(idx2, cA, sA);
        __m256d dx, dy; detail::residual4(d, precise, dx, dy);

        // s = dx*sA + dy*cA
    #if defined(__FMA__)
        __m256d s = _mm256_fmadd_pd(dx, sA, _mm256_mul_pd(dy, cA));
    #else
        __m256d s = _mm256_add_pd(_mm256_mul_pd(dx, sA), _mm256_mul_pd(dy, cA));
    #endif

        if (is_aligned_32(&sinOut[i])) _mm256_store_pd(&sinOut[i], s);
        else                            _mm256_storeu_pd(&sinOut[i], s);
    }

    // 2-lane SSE2 tail with AVX2 gathers
    for (; i + 1 < n; i += 2) {
        int     idx2[2];
        int64_t srem[2];

        for (int k = 0; k < 2; ++k) {
            int idx; uint64_t rem;
            btrig::phaseQ64_to_idx_rem(phase_q64[i + k], idx, rem);
            const uint64_t sign = rem >> 63;
            idx2[k] = (idx + (int)sign) & btrig::ANCHOR_MASK;
            srem[k] = (int64_t)rem;
        }

        const double K = (double)btrig::K_RESIDUAL;
        __m128d d = _mm_set_pd((double)srem[1] * K, (double)srem[0] * K);
        __m128d z = _mm_mul_pd(d, d);

        __m128i vi = _mm_loadu_si128(reinterpret_cast<const __m128i*>(idx2));
        __m128d cA = _mm_i32gather_pd(_anchorsCosD, vi, 8);
        __m128d sA = _mm_i32gather_pd(_anchorsSinD, vi, 8);

        __m128d dx, dy;
        if (precise) {
        #if defined(__FMA__)
            __m128d tc = _mm_set1_pd(-1.0/40320.0);
            tc = _mm_fmadd_pd(tc, z, _mm_set1_pd(1.0/720.0));
            tc = _mm_fmadd_pd(tc, z, _mm_set1_pd(-1.0/24.0));
            tc = _mm_fmadd_pd(tc, z, _mm_set1_pd(1.0/2.0));
            dx = _mm_fnmadd_pd(z, tc, _mm_set1_pd(1.0));

            __m128d ts = _mm_set1_pd(-1.0/5040.0);
            ts = _mm_fmadd_pd(ts, z, _mm_set1_pd(1.0/120.0));
            ts = _mm_fmadd_pd(ts, z, _mm_set1_pd(-1.0/6.0));
            ts = _mm_fmadd_pd(ts, z, _mm_set1_pd(1.0));
            dy = _mm_mul_pd(d, ts);
        #else
            __m128d tc = _mm_set1_pd(-1.0/40320.0);
            tc = _mm_add_pd(_mm_mul_pd(tc, z), _mm_set1_pd(1.0/720.0));
            tc = _mm_add_pd(_mm_mul_pd(tc, z), _mm_set1_pd(-1.0/24.0));
            tc = _mm_add_pd(tc, _mm_set1_pd(1.0/2.0));
            dx = _mm_sub_pd(_mm_set1_pd(1.0), _mm_mul_pd(z, tc));

            __m128d ts = _mm_set1_pd(-1.0/5040.0);
            ts = _mm_add_pd(_mm_mul_pd(ts, z), _mm_set1_pd(1.0/120.0));
            ts = _mm_add_pd(_mm_mul_pd(ts, z), _mm_set1_pd(-1.0/6.0));
            ts = _mm_add_pd(ts, _mm_set1_pd(1.0));
            dy = _mm_mul_pd(d, ts);
        #endif
        } else {
        #if defined(__FMA__)
            __m128d tc = _mm_set1_pd(-1.0/24.0);
            tc = _mm_fmadd_pd(tc, z, _mm_set1_pd(1.0/2.0));
            dx = _mm_fnmadd_pd(z, tc, _mm_set1_pd(1.0));

            __m128d ts = _mm_set1_pd(-1.0/6.0);
            ts = _mm_fmadd_pd(ts, z, _mm_set1_pd(1.0));
            dy = _mm_mul_pd(d, ts);
        #else
            __m128d tc = _mm_set1_pd(-1.0/24.0);
            tc = _mm_add_pd(_mm_mul_pd(tc, z), _mm_set1_pd(1.0/2.0));
            dx = _mm_sub_pd(_mm_set1_pd(1.0), _mm_mul_pd(z, tc));

            __m128d ts = _mm_set1_pd(-1.0/6.0);
            ts = _mm_add_pd(_mm_mul_pd(ts, z), _mm_set1_pd(1.0));
            dy = _mm_mul_pd(d, ts);
        #endif
        }

        // s = dx*sA + dy*cA
    #if defined(__FMA__)
        __m128d s = _mm_fmadd_pd(dx, sA, _mm_mul_pd(dy, cA));
    #else
        __m128d s = _mm_add_pd(_mm_mul_pd(dx, sA), _mm_mul_pd(dy, cA));
    #endif

        _mm_storeu_pd(&sinOut[i], s);
    }

    // Scalar tail
    for (; i < n; ++i) {
        long double c, s;
        btrig::sincos_q64(phase_q64[i], precise, c, s);
        sinOut[i] = (double)s;
    }
#else
    // No AVX2 at compile time -> scalar
    for (int i = 0; i < n; ++i) {
        long double c, s;
        btrig::sincos_q64(phase_q64[i], precise, c, s);
        sinOut[i] = (double)s;
    }
#endif
}


// -----------------------------------------------------------------------------
// Portable dispatcher: runtime AVX2 choose; otherwise scalar
// -----------------------------------------------------------------------------
static inline void sincos_q64(const uq64* phase_q64,
                              int n,
                              bool precise,
                              double* cosOut,
                              double* sinOut) {
    if (has_avx2()) {
        sincos_q64_avx2(phase_q64, n, precise, cosOut, sinOut);
        return;
    }
    for (int i = 0; i < n; ++i) {
        long double c, s;
        btrig::sincos_q64(phase_q64[i], precise, c, s);
        cosOut[i] = (double)c;
        sinOut[i] = (double)s;
    }
}

// Optional helpers that write long double outputs (SIMD path computes double then widens)
static inline void sincos_q64_ld(const uq64* phase_q64,
                                 int n,
                                 bool precise,
                                 long double* cosOut,
                                 long double* sinOut) {
    std::vector<double> cd(n), sd(n);
    sincos_q64(phase_q64, n, precise, cd.data(), sd.data());
    for (int i = 0; i < n; ++i) {
        cosOut[i] = (long double)cd[i];
        sinOut[i] = (long double)sd[i];
    }
}

// Simple wrappers that compute only cosine or only sine
static inline void cos_q64(const uq64* phase_q64, int n, bool precise, double* cosOut) {
    if (has_avx2()) { cos_q64_avx2(phase_q64, n, precise, cosOut); return; }
    for (int i=0; i<n; ++i) { long double c,s; btrig::sincos_q64(phase_q64[i], precise, c, s); cosOut[i]=(double)c; }
}
static inline void sin_q64(const uq64* phase_q64, int n, bool precise, double* sinOut) {
    if (has_avx2()) { sin_q64_avx2(phase_q64, n, precise, sinOut); return; }
    for (int i=0; i<n; ++i) { long double c,s; btrig::sincos_q64(phase_q64[i], precise, c, s); sinOut[i]=(double)s; }
}

}} // namespace btrig::simd

#endif // BTRIG_SIMD_H
