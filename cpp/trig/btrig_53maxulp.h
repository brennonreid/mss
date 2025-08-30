#pragma once
#ifndef BTRIG_H
#define BTRIG_H

#include <immintrin.h>
#include <cmath>    // for fma()

namespace btrig {

    // --------------------
    // Configuration
    // --------------------
    constexpr int TRIG_ANCHORS_BASE_POW = 5;
    constexpr int NUM_ANCHORS_QUADRANT = 1 << TRIG_ANCHORS_BASE_POW;
    constexpr int NUM_ANCHORS_FULL = NUM_ANCHORS_QUADRANT * 4;
    constexpr int ANCHOR_MASK = NUM_ANCHORS_FULL - 1;

    // --------------------
    // Core constants
    // --------------------
    constexpr double TAU = 6.283185307179586476925286766559; // 2*pi
    constexpr double QUADTAU = TAU / 4.0;
    constexpr double STEP = TAU / static_cast<double>(NUM_ANCHORS_FULL);
    constexpr double ONE_OVER_STEP = 1.0 / STEP;

    constexpr double DEG_STEP_QUADRANT = QUADTAU / static_cast<double>(NUM_ANCHORS_QUADRANT);
    constexpr double DEG_STEP_FULL = TAU / static_cast<double>(NUM_ANCHORS_FULL);

    // Precomputed inverse factorials for interpolation
    constexpr double INV_2 = 1.0 / 2.0;
    constexpr double INV_6 = 1.0 / 6.0;
    constexpr double INV_24 = 1.0 / 24.0;
    constexpr double INV_120 = 1.0 / 120.0;
    constexpr double INV_720 = 1.0 / 720.0;
    constexpr double INV_5040 = 1.0 / 5040.0;
    constexpr double INV_40320 = 1.0 / 40320.0;
    constexpr double INV_362880 = 2.755731922398589065255731922e-06; //  1/9!
    constexpr double INV_3628800 = 2.755731922398589065255731922e-07; // 1/10!

    // split-τ constants for hi/lo residual (Kahan-ish)
    constexpr double TAU_HIGH = 6.2831853071795862086997425149;
    constexpr double TAU_LOW = 2.449293598294706440262118471e-16;
    constexpr double STEP_HIGH = TAU_HIGH / static_cast<double>(NUM_ANCHORS_FULL);
    constexpr double STEP_LOW = TAU_LOW / static_cast<double>(NUM_ANCHORS_FULL);
    constexpr double INV_STEP = 1.0 / STEP_HIGH;

    constexpr double INV_TAU = 1.0 / TAU;

    // Anchor array: [angle][0=theta, 1=cos, 2=sin]
    static double fullCircleAnchors[NUM_ANCHORS_FULL][3];

    // --------------------
    // Internal helpers
    // --------------------
    static inline void applyQuadrantTransform(double inx, double iny, int quadrant,
        double& outx, double& outy) {
        switch (quadrant & 3) {
        case 0: outx = inx;  outy = iny;  break;
        case 1: outx = -iny; outy = inx;  break;
        case 2: outx = -inx; outy = -iny; break;
        default: outx = iny; outy = -inx; break;
        }
    }

    static inline double wrapTau(double a) {
        const double k = a * INV_TAU;
        long long q = (long long)k;
        if (a < 0.0 && (double)q != k) --q;
        return a - (double)q * TAU;
    }

    static inline double taylorCosFloat(double x, int depth) {
        double term = 1.0;
        double sum = 1.0;
        for (int i = 1; i <= depth; ++i) {
            term *= -x * x / ((2.0 * i - 1.0) * (2.0 * i));
            sum += term;
        }
        return sum;
    }

    static inline double taylorSinFloat(double x, int depth) {
        double term = x;
        double sum = x;
        for (int i = 1; i <= depth; ++i) {
            term *= -x * x / ((2.0 * i) * (2.0 * i + 1.0));
            sum += term;
        }
        return sum;
    }

    static inline void taylorSinCosFloat2(double angle, int depth,
        double& outSin, double& outCos) {
        angle = wrapTau(angle);
        int quadrant = static_cast<int>(angle / QUADTAU);
        double x = angle - quadrant * QUADTAU;

        double sinVal = taylorSinFloat(x, depth);
        double cosVal = taylorCosFloat(x, depth);

        applyQuadrantTransform(cosVal, sinVal, quadrant, outCos, outSin);
    }

    // round-to-nearest (ties-to-even) without lib calls
    static inline double rint_nearest(double x) {
        const double B = 6755399441055744.0;              // 2^52 + 2^51
        union { double d; unsigned long long u; } ux{ x }, ub{ B };
        ub.u |= (ux.u & 0x8000000000000000ULL);           // copy sign
        return (x + ub.d) - ub.d;
    }

    // Build anchors using reduction to r in [-pi/4, +pi/4]
    static void initFullCircleAnchors() {
        for (int i = 0; i < NUM_ANCHORS_FULL; ++i) {
            const double angle = i * DEG_STEP_FULL;

            // k = nearest multiple of (pi/2), r in [-pi/4, +pi/4]
            const double k_d = rint_nearest(angle / QUADTAU);   // QUADTAU = pi/2
            const long long k = (long long)k_d;
            const double r = angle - k_d * QUADTAU;             // |r| <= pi/4

            // small-argument series at r
            const double cos_r = taylorCosFloat(r, 21);
            const double sin_r = taylorSinFloat(r, 21);

            // rotate/flip back by k mod 4
            double c, s;
            switch (k & 3) {
            case 0:  c = cos_r;   s = sin_r;   break;
            case 1:  c = -sin_r;  s = cos_r;   break;
            case 2:  c = -cos_r;  s = -sin_r;  break;
            default: c = sin_r;  s = -cos_r;  break;
            }

            fullCircleAnchors[i][0] = angle;
            fullCircleAnchors[i][1] = c;
            fullCircleAnchors[i][2] = s;
        }
    }

    static struct _InitFullCircleAnchors {
        _InitFullCircleAnchors() { initFullCircleAnchors(); }
    } _initFullCircleAnchorsInstance;

    // --------------------
    // Custom trig funcs
    // --------------------
    static __forceinline void sin(double angle, bool precise, double& outSin) {
        // 1) index on the LUT grid (round-to-nearest, ties-to-even)
        const double t = angle * ONE_OVER_STEP;        // note: ONE_OVER_STEP, not 1/STEP_HIGH
        const double B = 6755399441055744.0;           // 2^52 + 2^51

        // sign-aware 2^52 trick (no lib call)
        union { double d; unsigned long long u; } tt{ t }, sb{ B };
        sb.u |= (tt.u & 0x8000000000000000ULL);        // copy sign of t onto B
        const double sB = sb.d;

        const double i_double = (t + sB) - sB;
        const int    idx = (int)i_double & ANCHOR_MASK;

        // 2) high-precision local residual: d = angle - i*(STEP_HIGH + STEP_LOW)
        const double step_mult_hi = i_double * STEP_HIGH;
        const double step_mult_lo = i_double * STEP_LOW;

        double d_hi = angle - step_mult_hi;
        double d_lo = -step_mult_lo;

        // single two-sum refine to keep tiny bits
        const double tmp = d_hi;
        d_hi = tmp + d_lo;
        d_lo = tmp - d_hi + d_lo;

        const double d = d_hi + d_lo;
        const double d2 = d * d;

        // 3) polynomial around anchor
        const double* a = fullCircleAnchors[idx];

        double dx, dy;
        if (precise) {
            dx = 1.0 - d2 * (INV_2 -
                d2 * (INV_24 -
                    d2 * (INV_720 -
                        d2 * INV_40320)));
            dy = d * (1.0 -
                d2 * (INV_6 -
                    d2 * (INV_120 -
                        d2 * INV_5040)));
        }
        else {
            dx = 1.0 - d2 * INV_2;
            dy = d;
        }

        // fused final combine: sin = dx*a2 + dy*a1
        outSin = fma(dy, a[1], dx * a[2]);
    }

    static __forceinline void cos(double angle, bool precise, double& outCos) {
        // 1) index on the LUT grid (round-to-nearest, ties-to-even)
        const double t = angle * ONE_OVER_STEP;        // note: ONE_OVER_STEP, not 1/STEP_HIGH
        const double B = 6755399441055744.0;           // 2^52 + 2^51

        // sign-aware 2^52 trick (no lib call)
        union { double d; unsigned long long u; } tt{ t }, sb{ B };
        sb.u |= (tt.u & 0x8000000000000000ULL);        // copy sign of t onto B
        const double sB = sb.d;

        const double i_double = (t + sB) - sB;
        const int    idx = (int)i_double & ANCHOR_MASK;

        // 2) high-precision local residual: d = angle - i*(STEP_HIGH + STEP_LOW)
        const double step_mult_hi = i_double * STEP_HIGH;
        const double step_mult_lo = i_double * STEP_LOW;

        double d_hi = angle - step_mult_hi;
        double d_lo = -step_mult_lo;

        // single two-sum refine to keep tiny bits
        const double tmp = d_hi;
        d_hi = tmp + d_lo;
        d_lo = tmp - d_hi + d_lo;

        const double d = d_hi + d_lo;
        const double d2 = d * d;

        // 3) polynomial around anchor
        const double* a = fullCircleAnchors[idx];

        double dx, dy;
        if (precise) {
            dx = 1.0 - d2 * (INV_2 -
                d2 * (INV_24 -
                    d2 * (INV_720 -
                        d2 * INV_40320)));
            dy = d * (1.0 -
                d2 * (INV_6 -
                    d2 * (INV_120 -
                        d2 * INV_5040)));
        }
        else {
            dx = 1.0 - d2 * INV_2;
            dy = d;
        }

        // fused final combine: cos = dx*a1 - dy*a2
        outCos = fma(dx, a[1], -dy * a[2]);
    }

    static __forceinline void sincos(double angle, bool precise,
        double& cosOut, double& sinOut)
    {
        // Use full step for indexing & thresholds
        const double STEP = STEP_HIGH + STEP_LOW;
        const double INV_STEP = 1.0 / STEP;

        // Nearest index (B-trick)
        const double t = fma(angle, INV_STEP, 0.0);
        const double B = 6755399441055744.0; // 2^52 + 2^51
        const double i_double = t + B - B;   // nearest integer in double
        const int64_t k0 = (int64_t)i_double; // safe for large angles
        int idx = (int)(k0 & ANCHOR_MASK);

        // Remainder for idx (hi/lo split)
        double d = fma(i_double, -STEP_HIGH, angle);
        d = fma(i_double, -STEP_LOW, d);

        // --- 3-way refine to the true nearest anchor ---
        {
            const double iL = i_double - 1.0;
            const double iR = i_double + 1.0;

            double dL = fma(iL, -STEP_HIGH, angle); dL = fma(iL, -STEP_LOW, dL);
            double dR = fma(iR, -STEP_HIGH, angle); dR = fma(iR, -STEP_LOW, dR);

            const double ad = fabs(d);
            const double adL = fabs(dL);
            const double adR = fabs(dR);

            if (adL < ad && adL <= adR) {
                idx = (int)((k0 - 1) & ANCHOR_MASK);
                d = dL;
            }
            else if (adR < ad && adR < adL) {
                idx = (int)((k0 + 1) & ANCHOR_MASK);
                d = dR;
            }
        }

        // Cardinal override (exact at multiples of quadrant size)
        const int Q = (NUM_ANCHORS_FULL >> 2);
        const bool on_cardinal = ((idx & (Q - 1)) == 0);

        const double* a = fullCircleAnchors[idx];
        double a1 = a[1], a2 = a[2];
        if (on_cardinal) {
            const int c = (idx >> TRIG_ANCHORS_BASE_POW) & 3; // 0,1,2,3 -> 0,90,180,270
            switch (c) {
            case 0:  a1 = 1.0; a2 = 0.0; break;
            case 1:  a1 = 0.0; a2 = 1.0; break;
            case 2:  a1 = -1.0; a2 = 0.0; break;
            default: a1 = 0.0; a2 = -1.0; break;
            }
        }

        // --- Taylor evaluation ---
        double dx, dy;
        if (precise) {
            const double d2 = d * d;

            // cos(d) up to d^10
            dx = fma(d2, -INV_3628800, INV_40320);
            dx = fma(d2, dx, -INV_720);
            dx = fma(d2, dx, INV_24);
            dx = fma(d2, dx, -INV_2);
            dx = fma(d2, dx, 1.0);

            // sin(d) up to d^9
            double p = fma(d2, INV_362880, -INV_5040);
            p = fma(d2, p, INV_120);
            p = fma(d2, p, -INV_6);
            p = fma(d2, p, 1.0);
            dy = d * p;
        }
        else {
            const double d2 = d * d;
            dx = fma(d2, -INV_2, 1.0);
            dy = d;
        }

        // --- Final rotation ---
        cosOut = fma(dx, a1, -dy * a2);
        sinOut = fma(dx, a2, dy * a1);
    }


    static __forceinline double atan2(double y, double x) {
        if (x == 0.0 && y == 0.0) return 0.0;

        const double ax = (y < 0.0 ? -y : y);
        const double ay = (x < 0.0 ? -x : x);

        const bool swap = (ax > ay);
        const double num = swap ? ay : ax;
        const double den = swap ? ax : ay;
        const double z = num / den;

        const double z2 = z * z;
        const double a = z * (0.9998660 +
            z2 * (-0.3302995 +
                z2 * (0.1801410 +
                    z2 * (-0.0851330 +
                        0.0208351 * z2))));

        const double angle = swap ? (1.5707963267948966 - a) : a;

        if (x >= 0.0) {
            return (y >= 0.0) ? angle : -angle;
        }
        else {
            return (y >= 0.0) ? (3.141592653589793 - angle) : (angle - 3.141592653589793);
        }
    }

    static_assert(sizeof(double) == 8, "fast_sqrt assumes 64-bit IEEE-754 double");

    static inline double fast_sqrt(double x) {
        if (x <= 0.0) return 0.0;
        union { double d; unsigned long long u; } v;
        v.d = x;
        v.u = 0x5fe6ec85e7de30daULL - (v.u >> 1);

        double y = v.d;
        const double xhalf = 0.5 * x;

        y = y * (1.5 - xhalf * y * y);
        y = y * (1.5 - xhalf * y * y);

        return x * y;
    }

    static inline double asin(double z) {
        const double HALFPI = TAU * 0.25;
        if (z >= 1.0)  return HALFPI;
        if (z <= -1.0) return -HALFPI;

        const double t = 1.0 - z * z;
        const double c = (t > 0.0) ? fast_sqrt(t) : 0.0;
        return atan2(z, c);
    }

    static inline double acos(double x) {
        if (x >= 1.0)  return 0.0;
        if (x <= -1.0) return TAU * 0.5;

        double t = 1.0 - x * x;
        double s = (t > 0.0) ? fast_sqrt(t) : 0.0;
        return atan2(s, x);
    }

    static __forceinline double tan(double angle) {
        double c, s;
        sincos(angle, /*precise=*/true, c, s);
        return s / c;
    }

    static __forceinline double atan(double y) {
        return atan2(y, 1.0);
    }

} // namespace btrig

#endif // BTRIG_H
