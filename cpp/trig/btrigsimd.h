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
        return _mm256_setr_pd(base[idx[0]], base[idx[1]], base[idx[2]], base[idx[3]]);
    }

    // Set to 0 to disable seam-bias completely
#ifndef BTRIG_EPS_SCALE
#define BTRIG_EPS_SCALE 0   // eps = STEP * 2^-BTRIG_EPS_SCALE
#endif

    static inline void sincos_batch_avx2(const double* angles, int n, bool precise,
        double* cosOut, double* sinOut) {
#if defined(__AVX2__)
        const __m256d ONE = _mm256_set1_pd(1.0);
        const __m256d ZERO = _mm256_set1_pd(0.0);
        const __m256d STEP_v = _mm256_set1_pd(btrig::STEP);
        const __m256d INV_STEP_v = _mm256_set1_pd(btrig::ONE_OVER_STEP);

        // use the same constants as scalar
        const __m256d INV_2_v = _mm256_set1_pd(btrig::INV_2);
        const __m256d INV_6_v = _mm256_set1_pd(btrig::INV_6);
        const __m256d INV_24_v = _mm256_set1_pd(btrig::INV_24);
        const __m256d INV_120_v = _mm256_set1_pd(btrig::INV_120);
        const __m256d INV_720_v = _mm256_set1_pd(btrig::INV_720);
        const __m256d INV_5040_v = _mm256_set1_pd(btrig::INV_5040);
        const __m256d INV_40320_v = _mm256_set1_pd(btrig::INV_40320);

        // modulo helper: r = floor(t) mod 2^k safely in double domain
        const double TWOPOW = double(btrig::ANCHOR_MASK + 1); // 2^k
        const __m256d TWOPOW_v = _mm256_set1_pd(TWOPOW);
        const __m256d INV_TWOPOW_v = _mm256_set1_pd(1.0 / TWOPOW);

        const __m128i MASK32 = _mm_set1_epi32(btrig::ANCHOR_MASK);

        int i = 0;
        for (; i + 3 < n; i += 4) {
            // t = angle * ONE_OVER_STEP
            __m256d ang = _mm256_loadu_pd(angles + i);
            __m256d t = _mm256_mul_pd(ang, INV_STEP_v);

            // tf = floor(t) - 1 (unconditional bias to previous cell, matching scalar)
            __m256d tf = _mm256_floor_pd(t);
            tf = _mm256_sub_pd(tf, ONE);

            // d = fma(-tf, STEP, angle)  (one rounding, same as scalar)
#if defined(__FMA__)
            __m256d d = _mm256_fmadd_pd(_mm256_sub_pd(ZERO, tf), STEP_v, ang);
#else
            __m256d d = _mm256_add_pd(ang, _mm256_mul_pd(_mm256_sub_pd(ZERO, tf), STEP_v));
#endif
            __m256d d2 = _mm256_mul_pd(d, d);

            // idx = (int)(tf mod 2^k)
            // q = floor(tf / 2^k); r = tf - q*2^k; idx = (int)r & mask
            __m256d q = _mm256_floor_pd(_mm256_mul_pd(tf, INV_TWOPOW_v));
            __m256d r = _mm256_sub_pd(tf, _mm256_mul_pd(q, TWOPOW_v));
            __m128i idx = _mm_and_si128(_mm256_cvttpd_epi32(r), MASK32);

            // micro-rotation polynomials (same order as scalar)
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
                // dx = 1 - d2/2; dy = d
                dx = _mm256_fnmadd_pd(d2, INV_2_v, ONE);
                dy = d;
            }

            // anchors (SoA, no gather)
            __m256d cosA = load4_pd_idx(A_COS, idx);
            __m256d sinA = load4_pd_idx(A_SIN, idx);

            // final rotation (same formulation as scalar)
#if defined(__FMA__)
            __m256d c = _mm256_fmsub_pd(dx, cosA, _mm256_mul_pd(dy, sinA)); // dx*c - dy*s
            __m256d s = _mm256_fmadd_pd(dx, sinA, _mm256_mul_pd(dy, cosA)); // dx*s + dy*c
#else
            __m256d c = _mm256_sub_pd(_mm256_mul_pd(dx, cosA), _mm256_mul_pd(dy, sinA));
            __m256d s = _mm256_add_pd(_mm256_mul_pd(dx, sinA), _mm256_mul_pd(dy, cosA));
#endif

            _mm256_storeu_pd(cosOut + i, c);
            _mm256_storeu_pd(sinOut + i, s);
        }

        // scalar tail (calls your exact scalar)
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
