#pragma once
#ifndef BTRIG_H
#define BTRIG_H

#include <immintrin.h>
#include <math.h> // for fma

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

    static void initFullCircleAnchors() {
        for (int i = 0; i < NUM_ANCHORS_FULL; ++i) {
            double angle = i * DEG_STEP_FULL;

            int quadrant = static_cast<int>(angle / QUADTAU);
            double x = fma(-(double)quadrant, QUADTAU, angle);

            double cosVal = taylorCosFloat(x, 22);
            double sinVal = taylorSinFloat(x, 22);

            // Used standard library cos/sin functions for higher accuracy
            // as an optimization trade-off, even though we'd prefer to use
            // our custom Taylor series-based methods for performance in the future.
            //double cosVal = std::cos(x);
            //double sinVal = std::sin(x);

            double cosFinal, sinFinal;
            applyQuadrantTransform(cosVal, sinVal, quadrant, cosFinal, sinFinal);

            fullCircleAnchors[i][0] = angle;
            fullCircleAnchors[i][1] = cosFinal;
            fullCircleAnchors[i][2] = sinFinal;
        }
    }

    static struct _InitFullCircleAnchors {
        _InitFullCircleAnchors() { initFullCircleAnchors(); }
    } _initFullCircleAnchorsInstance;

    // --------------------
    // Custom trig funcs
    // --------------------
    static __forceinline void sin(double angle, bool precise, double& outSin) {
        const double t = angle * ONE_OVER_STEP;
        int i = static_cast<int>(std::floor(t));

        const int idx = i & ANCHOR_MASK;

        const double d = std::fma(-(double)i, STEP, angle);
        const double d2 = d * d;                 // (mul is fine)
        const double* a = fullCircleAnchors[idx];

        double dx, dy;
        if (precise) {
            double p = std::fma(-d2, INV_40320, INV_720);
            p = std::fma(-d2, p, INV_24);
            p = std::fma(-d2, p, INV_2);
            dx = std::fma(-d2, p, 1.0);

            double q = std::fma(-d2, INV_5040, INV_120);
            q = std::fma(-d2, q, INV_6);
            dy = d * std::fma(-d2, q, 1.0);
        }
        else {
            dx = std::fma(-d2, INV_2, 1.0);
            dy = d;
        }

        // sin = dx*sin_anchor + dy*cos_anchor
        outSin = std::fma(dx, a[2], dy * a[1]);
    }

    static __forceinline void cos(double angle, bool precise, double& outCos) {
        const double t = angle * ONE_OVER_STEP;
        int i = static_cast<int>(std::floor(t));

        const int idx = i & ANCHOR_MASK;

        const double d = std::fma(-(double)i, STEP, angle);
        const double d2 = d * d;                 // (mul is fine)




        const double* a = fullCircleAnchors[idx];

        double dx, dy;
        if (precise) {
            double p = std::fma(-d2, INV_40320, INV_720);
            p = std::fma(-d2, p, INV_24);
            p = std::fma(-d2, p, INV_2);
            dx = std::fma(-d2, p, 1.0);

            double q = std::fma(-d2, INV_5040, INV_120);
            q = std::fma(-d2, q, INV_6);
            dy = d * std::fma(-d2, q, 1.0);
        }
        else {
            dx = std::fma(-d2, INV_2, 1.0);
            dy = d;
        }

        // cos = dx*cos_anchor - dy*sin_anchor
        outCos = std::fma(dx, a[1], -(dy * a[2]));
    }


    static __forceinline void sincos(double angle, bool precise,
        double& cosOut, double& sinOut) {
        // exact same i / idx logic as you have
        const double t = angle * ONE_OVER_STEP;
        int i = (int)t;
        i -= (t < 0.0);
        const int idx = i & ANCHOR_MASK;

        // compute residual in one rounding: d = angle - i*STEP
        const double d = std::fma(-(double)i, STEP, angle);
        const double d2 = d * d;                 // (mul is fine)


        const double* a = fullCircleAnchors[idx];

        double dx, dy;
        if (precise) {
            double p = std::fma(-d2, INV_40320, INV_720);
            p = std::fma(-d2, p, INV_24);
            p = std::fma(-d2, p, INV_2);
            dx = std::fma(-d2, p, 1.0);

            double q = std::fma(-d2, INV_5040, INV_120);
            q = std::fma(-d2, q, INV_6);
            dy = d * std::fma(-d2, q, 1.0);
        }
        else {
            dx = std::fma(-d2, INV_2, 1.0);
            dy = d;
        }

        cosOut = fma(dx, a[1], -(dy * a[2]));
        sinOut = fma(dx, a[2], (dy * a[1]));
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
