#pragma once
#include "btrig.h"

namespace btrig {

    // SoA anchor storage (one definition in a .cpp in production; inline here for brevity)
    alignas(64) double A_COS[NUM_ANCHORS_FULL];
    alignas(64) double A_SIN[NUM_ANCHORS_FULL];

    struct _InitSoAAnchors {
        _InitSoAAnchors() {
            for (int i = 0; i < NUM_ANCHORS_FULL; ++i) {
                A_COS[i] = fullCircleAnchors[i][1];
                A_SIN[i] = fullCircleAnchors[i][2];
            }
        }
    };
    static _InitSoAAnchors _initSoAAnchorsInstance;

    // helper: scalar-indexed "gather" for 4 doubles
    static inline __m256d load4_pd_idx(const double* base, __m128i idx32) {
        alignas(16) int idx[4];
        _mm_storeu_si128((__m128i*)idx, idx32);
        return _mm256_setr_pd(base[idx[0]], base[idx[1]], base[idx[2]], base[idx[3]]);
    }

    // -----------------------------------------------------------------------------
    // Precise-only, PHASE-native batch sincos (phase in turns, any real numbers OK)
    // -----------------------------------------------------------------------------
    static inline void sincos_batchp_avx2_precise(const double* phases, int n,
        double* cosOut, double* sinOut) {
#if defined(__AVX2__)
        const __m256d ONE = _mm256_set1_pd(1.0);
        const __m256d ZERO = _mm256_set1_pd(0.0);
        const __m256d N_v = _mm256_set1_pd((double)NUM_ANCHORS_FULL);
        const __m256d STEP_v = _mm256_set1_pd(btrig::STEP);

        // inverse factorials (Horner cos up to d^8, sin up to d^7)
        const __m256d INV_2_v = _mm256_set1_pd(btrig::INV_2);
        const __m256d INV_6_v = _mm256_set1_pd(btrig::INV_6);
        const __m256d INV_24_v = _mm256_set1_pd(btrig::INV_24);
        const __m256d INV_120_v = _mm256_set1_pd(btrig::INV_120);
        const __m256d INV_720_v = _mm256_set1_pd(btrig::INV_720);
        const __m256d INV_5040_v = _mm256_set1_pd(btrig::INV_5040);
        const __m256d INV_40320_v = _mm256_set1_pd(btrig::INV_40320);

        const __m128i MASK32 = _mm_set1_epi32(btrig::ANCHOR_MASK);

        int i = 0;
        for (; i + 3 < n; i += 4) {
            // Load 4 phases
            __m256d p = _mm256_loadu_pd(phases + i);

            // Normalize to [0,1): p = p - floor(p)
            __m256d pf = _mm256_floor_pd(p);
            p = _mm256_sub_pd(p, pf);

            // t = p * N; i = floor(t); r = t - i; idx = i & mask
            __m256d t = _mm256_mul_pd(p, N_v);
            __m256d tf = _mm256_floor_pd(t);

            __m128i i32 = _mm256_cvttpd_epi32(tf);              // safe since t >= 0
            __m256d i_pd = _mm256_cvtepi32_pd(i32);
            __m256d r = _mm256_sub_pd(t, i_pd);               // [0,1)

            __m128i idx = _mm_and_si128(i32, MASK32);

            // residual angle in radians: d = r * STEP
            __m256d d = _mm256_mul_pd(r, STEP_v);
            __m256d d2 = _mm256_mul_pd(d, d);

            // precise Horner: cos(d) ~ 1 - d^2/2 + d^4/24 - d^6/720 + d^8/40320
            __m256d z = _mm256_fnmadd_pd(d2, INV_40320_v, INV_720_v);
            z = _mm256_fnmadd_pd(d2, z, INV_24_v);
            __m256d in = _mm256_fnmadd_pd(d2, z, INV_2_v);
            __m256d dx = _mm256_fnmadd_pd(d2, in, ONE);

            // precise Horner: sin(d) ~ d*(1 - d^2/6 + d^4/120 - d^6/5040)
            __m256d u1 = _mm256_fnmadd_pd(d2, INV_5040_v, INV_120_v);
            __m256d u2 = _mm256_fnmadd_pd(d2, u1, INV_6_v);
            __m256d u3 = _mm256_fnmadd_pd(d2, u2, ONE);
            __m256d dy = _mm256_mul_pd(d, u3);

            // gather anchors
            __m256d ca = load4_pd_idx(A_COS, idx);
            __m256d sa = load4_pd_idx(A_SIN, idx);

            // final rotate:
            // cos = dx*ca - dy*sa;  sin = dx*sa + dy*ca
#if defined(__FMA__)
            __m256d c = _mm256_fmsub_pd(dx, ca, _mm256_mul_pd(dy, sa));
            __m256d s = _mm256_fmadd_pd(dx, sa, _mm256_mul_pd(dy, ca));
#else
            __m256d c = _mm256_sub_pd(_mm256_mul_pd(dx, ca), _mm256_mul_pd(dy, sa));
            __m256d s = _mm256_add_pd(_mm256_mul_pd(dx, sa), _mm256_mul_pd(dy, ca));
#endif

            _mm256_storeu_pd(cosOut + i, c);
            _mm256_storeu_pd(sinOut + i, s);
        }

        // tail
        for (; i < n; ++i) {
            double c, s;
            btrig::sincos(phases[i], /*precise=*/true, c, s);
            cosOut[i] = c; sinOut[i] = s;
        }
#else
        // Fallback: scalar precise
        for (int i = 0; i < n; ++i) {
            double c, s;
            btrig::sincos(phases[i], /*precise=*/true, c, s);
            cosOut[i] = c; sinOut[i] = s;
        }
#endif
    }

    // Optional convenience adapter (radians → phase) if you want it:
    static inline void sincos_batch_rad_precise(const double* radians, int n,
        double* cosOut, double* sinOut) {
        // Convert on the fly: phase = angle / TAU (no need to wrap; SIMD path wraps)
        // For best perf, preconvert to a phases array and call sincos_batchp_avx2_precise.
#if defined(__AVX2__)
        alignas(32) double tmp[4];
        int i = 0;
        const __m256d INV_TAU_v = _mm256_set1_pd(btrig::INV_TAU);
        for (; i + 3 < n; i += 4) {
            __m256d ang = _mm256_loadu_pd(radians + i);
            __m256d ph = _mm256_mul_pd(ang, INV_TAU_v);
            _mm256_store_pd(tmp, ph);
            sincos_batchp_avx2_precise(tmp, 4, cosOut + i, sinOut + i);
        }
        for (; i < n; ++i) {
            double p = radians[i] * btrig::INV_TAU;
            double c, s; btrig::sincos(p, /*precise=*/true, c, s);
            cosOut[i] = c; sinOut[i] = s;
        }
#else
        for (int i = 0; i < n; ++i) {
            double p = radians[i] * btrig::INV_TAU;
            double c, s; btrig::sincos(p, /*precise=*/true, c, s);
            cosOut[i] = c; sinOut[i] = s;
        }
#endif
    }

} // namespace btrig
