#pragma once

#include "btrig.h"


namespace btrig {

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

    // pack 4 indexed doubles (no gather; AMD-safe)
    static inline __m256d load4_pd_idx(const double* base, __m128i idx32) {
        alignas(16) int idx[4];
        _mm_storeu_si128((__m128i*)idx, idx32);
        return _mm256_setr_pd(base[idx[0]], base[idx[1]],
            base[idx[2]], base[idx[3]]);
    }

    static inline void sincos_batch_avx2(const double* angles, int n, bool precise,
        double* cosOut, double* sinOut) {
#if defined(__AVX2__)
        const __m256d ONE = _mm256_set1_pd(1.0);
        const __m256d HALF = _mm256_set1_pd(0.5);
        const __m256d STEP_v = _mm256_set1_pd(btrig::STEP);
        const __m256d INV_STEP_v = _mm256_set1_pd(btrig::ONE_OVER_STEP);

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
            // t = angle * ONE_OVER_STEP
            __m256d ang = _mm256_loadu_pd(angles + i);
            __m256d t = _mm256_mul_pd(ang, INV_STEP_v);

            // floor(t) and fractional offset
            __m256d tf = _mm256_floor_pd(t);
            __m256d d = _mm256_mul_pd(_mm256_sub_pd(t, tf), STEP_v);
            __m256d d2 = _mm256_mul_pd(d, d);

            // idx = (int)floor(t) & mask (32-bit indices)
            __m128i idx = _mm_and_si128(_mm256_cvttpd_epi32(tf), MASK32);

            // poly
            __m256d dx, dy;
            if (precise) {
                // dx = 1 - d2*(1/2 - d2*(1/24 - d2*(1/720 - d2*(1/40320))))
                __m256d z = _mm256_fnmadd_pd(d2, INV_40320_v, INV_720_v);
                z = _mm256_fnmadd_pd(d2, z, INV_24_v);
                __m256d in = _mm256_fnmadd_pd(d2, z, INV_2_v);
                dx = _mm256_fnmadd_pd(d2, in, ONE);
                // dy = d*(1 - d2*(1/6 - d2*(1/120 - d2*(1/5040))))
                __m256d u1 = _mm256_fnmadd_pd(d2, INV_5040_v, INV_120_v);
                __m256d u2 = _mm256_fnmadd_pd(d2, u1, INV_6_v);
                __m256d u3 = _mm256_fnmadd_pd(d2, u2, ONE);
                dy = _mm256_mul_pd(d, u3);
            }
            else {
                dx = _mm256_fnmadd_pd(d2, HALF, ONE); // 1 - 0.5*d2
                dy = d;
            }

            // anchors (no gather)
            __m256d cosA = load4_pd_idx(A_COS, idx);
            __m256d sinA = load4_pd_idx(A_SIN, idx);

            // rotate
#if defined(__FMA__)
            __m256d c = _mm256_fmsub_pd(dx, cosA, _mm256_mul_pd(dy, sinA));
            __m256d s = _mm256_fmadd_pd(dx, sinA, _mm256_mul_pd(dy, cosA));
#else
            __m256d c = _mm256_sub_pd(_mm256_mul_pd(dx, cosA), _mm256_mul_pd(dy, sinA));
            __m256d s = _mm256_add_pd(_mm256_mul_pd(dx, sinA), _mm256_mul_pd(dy, cosA));
#endif

            _mm256_storeu_pd(cosOut + i, c);
            _mm256_storeu_pd(sinOut + i, s);
        }

        // scalar tail
        for (; i < n; ++i) {
            double c, s;
            btrig::sincos(angles[i], precise, c, s);
            cosOut[i] = c;
            sinOut[i] = s;
        }
#else
        // scalar fallback
        for (int i = 0; i < n; ++i) {
            double c, s;
            btrig::sincos(angles[i], precise, c, s);
            cosOut[i] = c;
            sinOut[i] = s;
        }
#endif
    }

} // namespace btrig
