#pragma once
#ifndef BTRIG_H
#define BTRIG_H

#include <math.h> // for fma
#include <immintrin.h>

namespace btrig {

    // --------------------
    // Configuration
    // --------------------
    constexpr int TRIG_ANCHORS_BASE_POW = 6;
    constexpr int NUM_ANCHORS_QUADRANT = 1 << TRIG_ANCHORS_BASE_POW;
    constexpr int NUM_ANCHORS_FULL = NUM_ANCHORS_QUADRANT * 4;
    constexpr int ANCHOR_MASK = NUM_ANCHORS_FULL - 1;

    // --------------------
    // Core constants
    // --------------------
    constexpr double TAU = 6.283185307179586476925286766559;
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
    constexpr double INV_362880 = 2.755731922398589065255731922e-06;
    constexpr double INV_3628800 = 2.755731922398589065255731922e-07;

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
        return fma(-(double)q, TAU, a);
    }

    static inline double taylorCosFloat(double x, int depth) {
        double term = 1.0;
        double sum = 1.0;
        for (int i = 1; i <= depth; ++i) {
            term *= -x * x / ((2.0 * i - 1.0) * (2.0 * i));
            sum = fma(1.0, term, sum); // sum += term (via fma to fuse mults inside term if compiled aggressively)
        }
        return sum;
    }

    static inline double taylorSinFloat(double x, int depth) {
        double term = x;
        double sum = x;
        for (int i = 1; i <= depth; ++i) {
            term *= -x * x / ((2.0 * i) * (2.0 * i + 1.0));
            sum = fma(1.0, term, sum); // sum += term
        }
        return sum;
    }

    static inline void taylorSinCosFloat2(double angle, int depth,
        double& outSin, double& outCos) {
        angle = wrapTau(angle);
        int quadrant = static_cast<int>(angle / QUADTAU);
        double x = fma(-(double)quadrant, QUADTAU, angle);

        double sinVal = taylorSinFloat(x, depth);
        double cosVal = taylorCosFloat(x, depth);

        applyQuadrantTransform(cosVal, sinVal, quadrant, outCos, outSin);
    }

    // Function using long double for higher precision during initialization
    static inline long double taylorCosLongDouble(double x, int depth) {
        long double term = 1.0L;  // long double for higher precision
        long double sum = 1.0L;
        for (int i = 1; i <= depth; ++i) {
            term *= fma(-x * x, 1.0L / ((2.0L * i - 1.0L) * (2.0L * i)), 0.0L);
            // using long double
            sum = fma(1.0L, term, sum); // sum += term (via fma to fuse multiplications inside term)
        }
        return sum;
    }

    static inline long double taylorSinLongDouble(double x, int depth) {
        long double term = x, sum = x;
        const long double x2 = (long double)x * (long double)x;  // optional: hoist
        for (int i = 1; i <= depth; ++i) {
            term *= fma(-x2, 1.0L / ((2.0L * i) * (2.0L * i + 1.0L)), 0.0L);
            sum = fma(1.0L, term, sum);
        }
        return sum;
    }

    // Evaluate cos(x) with Horner in long double using even Taylor terms up to max_even (e.g. 22)
    static inline long double cos_horner_ld(double x, int max_even) {
        // Highest even power actually used
        const int twoN = (max_even & ~1);         // floor to even
        const int N = twoN / 2;

        const long double z = (long double)x * (long double)x; // x^2

        // c_N = (-1)^N / (2N)!
        long double cN = 1.0L;
        for (int m = 2; m <= 2 * N; ++m) cN /= (long double)m;
        if (N & 1) cN = -cN;

        long double acc = cN;
        // Horner down: c_k = -c_{k+1} * (2k+2)(2k+1)
        for (int k = N - 1; k >= 1; --k) {
            cN = -cN * (long double)(2 * (k + 1)) * (long double)(2 * (k + 1) - 1);
            acc = fma(acc, z, cN);
        }
        // Final "+ c0" with c0 = 1
        acc = fma(acc, z, 1.0L);
        return acc;
    }

    // Evaluate sin(x) with Horner in long double using odd Taylor terms up to max_odd (e.g. 21 or 23)
    static inline long double sin_horner_ld(double x, int max_odd) {
        // Highest odd power actually used
        const int twoNp1 = (max_odd | 1);         // ceil to odd
        const int N = (twoNp1 - 1) / 2;

        const long double xl = (long double)x;
        const long double z = xl * xl;           // x^2

        // s_N = (-1)^N / (2N+1)!
        long double sN = 1.0L;
        for (int m = 2; m <= 2 * N + 1; ++m) sN /= (long double)m;
        if (N & 1) sN = -sN;

        long double acc = sN;
        // Horner down: s_k = -s_{k+1} * (2k+3)(2k+2)
        for (int k = N - 1; k >= 0; --k) {
            sN = -sN * (long double)(2 * (k + 1) + 1) * (long double)(2 * (k + 1));
            acc = fma(acc, z, sN);
        }
        // sin(x) = x * P(z)
        return fma(xl, acc, 0.0L);
    }

    // Convenience: init-time sin/cos via Horner (keep your quadrant handling)
    static inline void hornerSinCosInit(double angle, int depth, double& outSin, double& outCos) {
        // Your wrap + quadrant split
        angle = wrapTau(angle);
        int quadrant = (int)(angle / QUADTAU);
        double x = fma(-(double)quadrant, QUADTAU, angle);

        long double c = cos_horner_ld(x, depth);         // even terms up to 'depth'
        long double s = sin_horner_ld(x, depth - 1);     // odd terms up to 'depth-1' (typical pairing)

        //long double c = std::cos(x);
        //long double s = std::sin(x);

        double cosVal = (double)c;
        double sinVal = (double)s;
        applyQuadrantTransform(cosVal, sinVal, quadrant, outCos, outSin);
    }

    // …and in your LUT builder:
    static void initFullCircleAnchors() {
        for (int i = 0; i < NUM_ANCHORS_FULL; ++i) {
            double angle = i * DEG_STEP_FULL;
            double s, c;
            
            hornerSinCosInit(angle, /*depth=*/22, s, c);
            //s = std::sin(angle);
            ///c = std::cos(angle);
            
            fullCircleAnchors[i][0] = angle;
            fullCircleAnchors[i][1] = c;
            fullCircleAnchors[i][2] = s;
        }
    }

    /*
    // Full Circle Table Initialization
    static void initFullCircleAnchors() {
        for (int i = 0; i < NUM_ANCHORS_FULL; ++i) {
            double angle = i * DEG_STEP_FULL;

            int quadrant = static_cast<int>(angle / QUADTAU);
            double x = fma(-(double)quadrant, QUADTAU, angle);

            // Using high precision during table initialization
            long double cosVal = taylorCosLongDouble(x, 12);  // Higher precision initialization
            long double sinVal = taylorSinLongDouble(x, 12);  // Higher precision initialization

            // Apply the quadrant transform and store the result
            double cosFinal, sinFinal;
            applyQuadrantTransform(static_cast<double>(cosVal), static_cast<double>(sinVal), quadrant, cosFinal, sinFinal);

            // Store the result back into the LUT as double
            fullCircleAnchors[i][0] = angle;
            fullCircleAnchors[i][1] = cosFinal;
            fullCircleAnchors[i][2] = sinFinal;
        }
    }
    */

    static struct _InitFullCircleAnchors {
        _InitFullCircleAnchors() { initFullCircleAnchors(); }
    } _initFullCircleAnchorsInstance;

    // --------------------
    // Custom trig funcs
    // --------------------
    static __forceinline void sin(double angle, bool precise, double& outSin) {
        const double t = angle * ONE_OVER_STEP;
        int i = static_cast<int>(floor(t));

        // seam-cancel pivot: previous cell
        const int idx = (i - 1) & ANCHOR_MASK;

        // residual from that pivot: angle - (i-1)*STEP  => in [STEP, 2*STEP)
        const double d = fma(-(double)(i - 1), STEP, angle);
        const double d2 = d * d;
        const double* a = fullCircleAnchors[idx];

        double dx, dy;
        if (precise) {
            double p = fma(-d2, INV_40320, INV_720);
            p = fma(-d2, p, INV_24);
            p = fma(-d2, p, INV_2);
            dx = fma(-d2, p, 1.0);

            double q = fma(-d2, INV_5040, INV_120);
            q = fma(-d2, q, INV_6);
            dy = d * fma(-d2, q, 1.0);
        }
        else {
            dx = fma(-d2, INV_2, 1.0);
            dy = d;
        }

        // sin = dx*sin_anchor + dy*cos_anchor
        outSin = fma(dx, a[2], dy * a[1]);
    }

    static __forceinline void cos(double angle, bool precise, double& outCos) {
        const double t = angle * ONE_OVER_STEP;
        int i = static_cast<int>(floor(t));

        // seam-cancel pivot: previous cell
        const int idx = (i - 1) & ANCHOR_MASK;

        // residual from that pivot
        const double d = fma(-(double)(i - 1), STEP, angle);
        const double d2 = d * d;
        const double* a = fullCircleAnchors[idx];

        double dx, dy;
        if (precise) {
            double p = fma(-d2, INV_40320, INV_720);
            p = fma(-d2, p, INV_24);
            p = fma(-d2, p, INV_2);
            dx = fma(-d2, p, 1.0);

            double q = fma(-d2, INV_5040, INV_120);
            q = fma(-d2, q, INV_6);
            dy = d * fma(-d2, q, 1.0);
        }
        else {
            dx = fma(-d2, INV_2, 1.0);
            dy = d;
        }

        // cos = dx*cos_anchor - dy*sin_anchor
        outCos = fma(dx, a[1], -(dy * a[2]));
    }

    static __forceinline void sincos(double angle, bool precise,
        double& cosOut, double& sinOut)
    {
        const double t = angle * ONE_OVER_STEP;
        int i = static_cast<int>(floor(t));

        // seam-cancel pivot: previous cell
        const int idx = (i - 1) & ANCHOR_MASK;

        // residual from that pivot
        const double d = fma(-(double)(i - 1), STEP, angle);
        const double d2 = d * d;
        const double* a = fullCircleAnchors[idx];

        double dx, dy;
        if (precise) {
            double p = fma(-d2, INV_40320, INV_720);
            p = fma(-d2, p, INV_24);
            p = fma(-d2, p, INV_2);
            dx = fma(-d2, p, 1.0);

            double q = fma(-d2, INV_5040, INV_120);
            q = fma(-d2, q, INV_6);
            dy = d * fma(-d2, q, 1.0);
        }
        else {
            dx = fma(-d2, INV_2, 1.0);
            dy = d;
        }

        // rotate about biased anchor
        cosOut = fma(dx, a[1], -(dy * a[2]));
        sinOut = fma(dx, a[2], dy * a[1]);
    }





    static __forceinline double tan(double angle) {
        double c, s;
        sincos(angle, /*precise=*/true, c, s);
        return s / c;
    }



    /*
    static_assert(sizeof(double) == 8, "fast_sqrt assumes 64-bit IEEE-754 double");

    static inline double fast_sqrt(double x) {
        if (x <= 0.0) return 0.0;
        union { double d; unsigned long long u; } v;
        v.d = x;
        v.u = 0x5fe6ec85e7de30daULL - (v.u >> 1);

        double y = v.d;
        const double xhalf = 0.5 * x;

        y = y * (1.5 - fma(xhalf, y * y, 0.0)); // y *= (1.5 - xhalf*y*y)
        y = y * (1.5 - fma(xhalf, y * y, 0.0));

        return x * y;
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
        double p = fma(0.0208351, z2, -0.0851330);
        p = fma(p, z2, 0.1801410);
        p = fma(p, z2, -0.3302995);
        p = fma(p, z2, 0.9998660);
        const double a = z * p;

        const double angle = swap ? (1.5707963267948966 - a) : a;

        if (x >= 0.0) {
            return (y >= 0.0) ? angle : -angle;
        }
        else {
            return (y >= 0.0) ? (3.141592653589793 - angle) : (angle - 3.141592653589793);
        }
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

    static __forceinline double atan(double y) {
        return atan2(y, 1.0);
    }
    */
} // namespace btrig

#endif // BTRIG_H
