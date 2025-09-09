#pragma once
#ifndef BTRIG_H
#define BTRIG_H

#include <cmath>
#include <immintrin.h>


namespace btrig {

    // -------------------------------------------------------------------------
    // Configuration
    // -------------------------------------------------------------------------
    // Phase-driven trig with a uniform anchor table over one full turn.
    // TRIG_ANCHORS_FULL_POW controls granularity: 2^pow anchors per full turn.
    // Total anchors across the full turn = 2^pow.
    constexpr int TRIG_ANCHORS_FULL_POW = 8;

    // Calculate the total number of anchors in one step using the power.
    constexpr int NUM_ANCHORS_FULL = 1 << TRIG_ANCHORS_FULL_POW;

    // The mask remains the same, used for bitwise lookups.
    constexpr int ANCHOR_MASK = NUM_ANCHORS_FULL - 1;

    // -------------------------------------------------------------------------
    // Core constants
    // -------------------------------------------------------------------------
    // Note: The PUBLIC API uses phase ("turns"): 1.0 == full rotation (TAU radians).
    // Some internals (anchor init, inverse trig) temporarily use radians, then
    // convert back to phase; that conversion is strictly internal.
    constexpr double TAU = 6.283185307179586476925286766559;
    constexpr double PI = TAU / 2.0;
    constexpr double PI_2 = TAU / 4.0;
    constexpr double PI_4 = TAU / 8.0;

    constexpr double QUADTAU = TAU / 4.0;

    // STEP is the angle in radians between adjacent anchors.
    // ONE_OVER_STEP is provided for optional callers who want that ratio.
    constexpr double STEP = TAU / static_cast<double>(NUM_ANCHORS_FULL);
    constexpr double ONE_OVER_STEP = 1.0 / STEP;

    // Precomputed inverse factorials used by small-angle Taylor/Horner evaluation.
    constexpr double INV_2 = 1.0 / 2.0;
    constexpr double INV_6 = 1.0 / 6.0;
    constexpr double INV_24 = 1.0 / 24.0;
    constexpr double INV_120 = 1.0 / 120.0;
    constexpr double INV_720 = 1.0 / 720.0;
    constexpr double INV_5040 = 1.0 / 5040.0;
    constexpr double INV_40320 = 1.0 / 40320.0;
    //constexpr double INV_362880 = 2.755731922398589065255731922e-06;
    //constexpr double INV_3628800 = 2.755731922398589065255731922e-07;

    // Phase <-> radians conversion factor for internal use.
    constexpr double INV_TAU = 1.0 / TAU;

    // Anchor array:
    //   per anchor i:
    //     [0] = phase in turns (i / NUM_ANCHORS_FULL)
    //     [1] = cos(angle) at that anchor
    //     [2] = sin(angle) at that anchor
    //
    // At call-time, we pick a pivot anchor and rotate it by a small residual.
    static double fullCircleAnchors[NUM_ANCHORS_FULL][3];

    // -------------------------------------------------------------------------
    // Internal helpers
    // -------------------------------------------------------------------------
    // Quadrant transform helper kept for angle-based experiments/init paths.
    static inline void applyQuadrantTransform(double inx, double iny, int quadrant,
        double& outx, double& outy) {
        switch (quadrant & 3) {
        case 0: outx = inx;  outy = iny;  break;
        case 1: outx = -iny; outy = inx;  break;
        case 2: outx = -inx; outy = -iny; break;
        default: outx = iny; outy = -inx; break;
        }
    }

    // Legacy/radian helper used by internal angle-based functions.
    // Converts an angle (radians) into [-TAU, TAU) repeatedly until within [0, TAU),
    // implemented as a turn-scaled floor and a fused correction.
    static inline double wrapTau(double a) {
        // Scale into turns
        double k = a * INV_TAU;
        // Round toward -∞ without int conversion
        double q = floor(k);
        return fma(-q, TAU, a);
    }

    // Phase wrapper: returns fractional part in [0,1) for any real input.
    // Provided for optional external use; sin/cos/sincos already self-wrap.
    static __forceinline double wrapPhase(double p) {
        double ip, frac = std::modf(p, &ip);     // frac in (-1,1)
        return (frac >= 0.0) ? frac : (frac + 1.0);
    }


    // -------------------------------------------------------------------------
    // Angle-based Taylor helpers (kept for experiments / anchor init only).
    // These accept radians by design and are not called in the phase API paths.
    // -------------------------------------------------------------------------

    static inline double taylorCosFloat(double x, int depth) {
        double term = 1.0;
        double sum = 1.0;
        for (int i = 1; i <= depth; ++i) {
            term *= -x * x / ((2.0 * i - 1.0) * (2.0 * i));
            sum = fma(1.0, term, sum); // accumulate with FMA
        }
        return sum;
    }

    static inline double taylorSinFloat(double x, int depth) {
        double term = x;
        double sum = x;
        for (int i = 1; i <= depth; ++i) {
            term *= -x * x / ((2.0 * i) * (2.0 * i + 1.0));
            sum = fma(1.0, term, sum); // accumulate with FMA
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

    // Long-double Horner evaluators for experiments/init (radians in, not used by phase API).
    static inline long double taylorCosLongDouble(double x, int depth) {
        long double term = 1.0L;
        long double sum = 1.0L;
        for (int i = 1; i <= depth; ++i) {
            term *= fma(-x * x, 1.0L / ((2.0L * i - 1.0L) * (2.0L * i)), 0.0L);
            sum = fma(1.0L, term, sum);
        }
        return sum;
    }

    static inline long double taylorSinLongDouble(double x, int depth) {
        long double term = x, sum = x;
        const long double x2 = (long double)x * (long double)x;
        for (int i = 1; i <= depth; ++i) {
            term *= fma(-x2, 1.0L / ((2.0L * i) * (2.0L * i + 1.0L)), 0.0L);
            sum = fma(1.0L, term, sum);
        }
        return sum;
    }

    // Horner (even/odd) cos/sin in long double, for optional init/validation use.
    static inline long double cos_horner_ld(double x, int max_even) {
        const int twoN = (max_even & ~1);
        const int N = twoN / 2;

        const long double z = (long double)x * (long double)x; // x^2

        long double cN = 1.0L;
        for (int m = 2; m <= 2 * N; ++m) cN /= (long double)m;
        if (N & 1) cN = -cN;

        long double acc = cN;
        for (int k = N - 1; k >= 1; --k) {
            cN = -cN * (long double)(2 * (k + 1)) * (long double)(2 * (k + 1) - 1);
            acc = fma(acc, z, cN);
        }
        acc = fma(acc, z, 1.0L);
        return acc;
    }

    static inline long double sin_horner_ld(double x, int max_odd) {
        const int twoNp1 = (max_odd | 1);
        const int N = (twoNp1 - 1) / 2;

        const long double xl = (long double)x;
        const long double z = xl * xl;

        long double sN = 1.0L;
        for (int m = 2; m <= 2 * N + 1; ++m) sN /= (long double)m;
        if (N & 1) sN = -sN;

        long double acc = sN;
        for (int k = N - 1; k >= 0; --k) {
            sN = -sN * (long double)(2 * (k + 1) + 1) * (long double)(2 * (k + 1));
            acc = fma(acc, z, sN);
        }
        return fma(xl, acc, 0.0L);
    }

    // Convenience: init-time sin/cos via Horner (angle/radian path).
    // Kept for experiments; not used by the phase API at runtime.
    static inline void hornerSinCosInit(double angle, int depth, double& outSin, double& outCos) {
        angle = wrapTau(angle);
        int quadrant = (int)(angle / QUADTAU);
        double x = fma(-(double)quadrant, QUADTAU, angle);

        long double c = cos_horner_ld(x, depth);
        long double s = sin_horner_ld(x, depth - 1);

        //long double c = std::cos(x);
        //long double s = std::sin(x);

        double cosVal = (double)c;
        double sinVal = (double)s;
        applyQuadrantTransform(cosVal, sinVal, quadrant, outCos, outSin);
    }

    // -------------------------------------------------------------------------
    // Anchor initialization
    // -------------------------------------------------------------------------
    // Two variants:
    //  - initFullCircleAnchors(): directly computes anchors via std::sin/cos on angle=phase*TAU.
    //  - initFullCircleAnchors2(): iteratively rotates by STEP using a fixed rotation; kept as an
    //    alternative for drift experiments (currently not selected).
    //
    // In both cases we store:
    //   [0]=phase (turns), [1]=cos(anchor_angle), [2]=sin(anchor_angle)
    // and we pin exact cardinals by construction for seam robustness.
    static void initFullCircleAnchors3() {
        // STEP is in radians between anchors; here angle = phase*TAU only to compute cos/sin at init.
        const double step = STEP;

        for (int i = 0; i < NUM_ANCHORS_FULL; ++i) {
            double phase = (double)i / NUM_ANCHORS_FULL;  // exact fraction of a turn
            double angle = phase * TAU;                   // internal: angle used once here
            fullCircleAnchors[i][0] = phase;              // store phase
            fullCircleAnchors[i][1] = std::cos(angle);
            fullCircleAnchors[i][2] = std::sin(angle);
        }

        // Force exact cardinal points
        constexpr int Q = NUM_ANCHORS_FULL / 4;

        fullCircleAnchors[0][1] = 1.0;   fullCircleAnchors[0][2] = 0.0;
        fullCircleAnchors[1 * Q][1] = 0.0;   fullCircleAnchors[1 * Q][2] = 1.0;
        fullCircleAnchors[2 * Q][1] = -1.0;  fullCircleAnchors[2 * Q][2] = 0.0;
        fullCircleAnchors[3 * Q][1] = 0.0;   fullCircleAnchors[3 * Q][2] = -1.0;
        //fullCircleAnchors[NUM_ANCHORS_FULL][1] = 1.0; // wrap back to start
        //fullCircleAnchors[NUM_ANCHORS_FULL][2] = 0.0;
    }


    // --------------------
// Anchor initialization (Horner-only; no std::sin/cos)
// --------------------
    static void initFullCircleAnchors() {
        // Choose an even depth for cos; sin uses (depth-1).
        // 22/21 matches your earlier high-accuracy setup.
        constexpr int HORNER_DEPTH = 22;

        for (int i = 0; i < NUM_ANCHORS_FULL; ++i) {
            const double phase = (double)i / NUM_ANCHORS_FULL;  // exact fraction of a turn
            const double angle = phase * TAU;                   // radians for Horner eval

            double s, c;
            hornerSinCosInit(angle, HORNER_DEPTH, s, c);

            fullCircleAnchors[i][0] = phase; // store phase (turns)
            fullCircleAnchors[i][1] = c;     // cos(angle)
            fullCircleAnchors[i][2] = s;     // sin(angle)
        }

        // Force exact cardinal points for perfect seams.
        constexpr int Q = NUM_ANCHORS_FULL / 4;
        fullCircleAnchors[0][1] = 1.0;  fullCircleAnchors[0][2] = 0.0;  // 0
        fullCircleAnchors[1 * Q][1] = 0.0;  fullCircleAnchors[1 * Q][2] = 1.0;  // π/2
        fullCircleAnchors[2 * Q][1] = -1.0; fullCircleAnchors[2 * Q][2] = 0.0;  // π
        fullCircleAnchors[3 * Q][1] = 0.0;  fullCircleAnchors[3 * Q][2] = -1.0; // 3π/2
    }


    // Static initializer chooses the direct std::sin/cos anchor init above.
    static struct _InitFullCircleAnchors {
        _InitFullCircleAnchors() { initFullCircleAnchors(); }
    } _initFullCircleAnchorsInstance;

    // -------------------------------------------------------------------------
    // Phase-based trig API
    // -------------------------------------------------------------------------
    //
    // PHASE semantics:
    //   - Inputs are in "turns": any real number is accepted.
    //     The functions normalize phase internally to [0,1) using floor-based wrap.
    //   - 1.0 turn == TAU radians (full rotation).
    //   - The algorithm selects a pivot anchor (i-1), computes a small residual
    //     angle d in radians (in [STEP, 2*STEP)), evaluates short Taylor/Horner
    //     series for cos(d)/sin(d), and rotates the pivot anchor using FMA.
    //
    // precise=false  : very fast small-angle approx (cos ~ 1 - d^2/2, sin ~ d)
    // precise=true   : higher-order Horner (cos to 8th, sin to 7th) for high accuracy
    //
    // All functions are robust at seam/cardinal/diagonal boundaries by construction.
    static __forceinline void sincos(double phase, bool precise,
        double& cosOut, double& sinOut)
    {
        double p = phase - std::floor(phase); // normalize to [0,1)

        // Anchor index and residual
        const double t = p * NUM_ANCHORS_FULL;     // [0, N)
        const int    i = (int)t;                   // floor(t)
        const int    idx = (i - 1) & ANCHOR_MASK;  // pivot anchor (i-1)
        const double* a = fullCircleAnchors[idx];  // [phase, cos, sin]

        // residual angle d (radians) relative to pivot
        const double r = t - (double)i;           // [0,1)
        const double d = std::fma(r, STEP, STEP); // (r + 1.0) * STEP
        const double d2 = d * d;

        double dx, dy;

        if (precise) {
            // Cos(d) ≈ 1 - d^2/2! + d^4/4! - d^6/6! + d^8/8!  (Horner with FMA)
            double poly = std::fma(-INV_40320, d2, INV_720);
            poly = std::fma(-d2, poly, INV_24);
            poly = std::fma(-d2, poly, INV_2);
            dx = std::fma(-d2, poly, 1.0);

            // Sin(d) ≈ d * [1 - d^2/3! + d^4/5! - d^6/7!]
            double q = std::fma(-INV_5040, d2, INV_120);
            q = std::fma(-d2, q, INV_6);
            dy = d * std::fma(-d2, q, 1.0);
        }
        else {
            // Fast small-angle approx
            dx = std::fma(-0.5, d2, 1.0); // 1 - 0.5*d^2
            dy = d;
        }

        // Rotate pivot anchor by residual (FMA to reduce rounding)
        sinOut = std::fma(dx, a[2], dy * a[1]); // sa*dx + ca*dy
        cosOut = std::fma(dx, a[1], -dy * a[2]); // ca*dx - sa*dy
    }


    // Singleton fast paths (phase API)
    static __forceinline void sin(double phase, bool precise, double& sinOut) {
        double p = phase - std::floor(phase); // [0,1)

        const double t = p * NUM_ANCHORS_FULL;
        const int    i = (int)t;
        const int    idx = (i - 1) & ANCHOR_MASK;
        const double* a = fullCircleAnchors[idx];

        const double r = t - (double)i;
        const double d = std::fma(r, STEP, STEP);
        const double d2 = d * d;

        double dx, dy;

        if (precise) {
            double poly = std::fma(-INV_40320, d2, INV_720);
            poly = std::fma(-d2, poly, INV_24);
            poly = std::fma(-d2, poly, INV_2);
            dx = std::fma(-d2, poly, 1.0);

            double q = std::fma(-INV_5040, d2, INV_120);
            q = std::fma(-d2, q, INV_6);
            dy = d * std::fma(-d2, q, 1.0);
        }
        else {
            dx = std::fma(-0.5, d2, 1.0);  // 1 - 0.5*d^2
            dy = d;
        }

        sinOut = std::fma(dx, a[2], dy * a[1]);  // sa*dx + ca*dy
    }

    static __forceinline void cos(double phase, bool precise, double& cosOut) {
        double p = phase - std::floor(phase); // [0,1)

        const double t = p * NUM_ANCHORS_FULL;
        const int    i = (int)t;
        const int    idx = (i - 1) & ANCHOR_MASK;
        const double* a = fullCircleAnchors[idx];

        const double r = t - (double)i;
        const double d = std::fma(r, STEP, STEP);
        const double d2 = d * d;

        double dx, dy;

        if (precise) {
            double poly = std::fma(-INV_40320, d2, INV_720);
            poly = std::fma(-d2, poly, INV_24);
            poly = std::fma(-d2, poly, INV_2);
            dx = std::fma(-d2, poly, 1.0);

            double q = std::fma(-INV_5040, d2, INV_120);
            q = std::fma(-d2, q, INV_6);
            dy = d * std::fma(-d2, q, 1.0);
        }
        else {
            dx = std::fma(-0.5, d2, 1.0);  // 1 - 0.5*d^2
            dy = d;
        }

        cosOut = std::fma(dx, a[1], -dy * a[2]);  // ca*dx - sa*dy
    }

    // tan from phase: uses precise sincos for stability near cos≈0
    static __forceinline double tan(double phase) {
        double c, s;
        sincos(phase, /*precise=*/true, c, s);
        return s / c;
    }

    // -------------------------------------------------------------------------
    // Value-based functions (kept and used internally)
    // -------------------------------------------------------------------------
    static_assert(sizeof(double) == 8, "fast_sqrt assumes 64-bit IEEE-754 double");

    // Fast sqrt via reciprocal-sqrt Newton steps, using FMA inside refinement.
    static inline double fast_sqrt(double x) {
        if (x <= 0.0) return 0.0;
        union { double d; unsigned long long u; } v;
        v.d = x;
        v.u = 0x5fe6ec85e7de30daULL - (v.u >> 1);

        double y = v.d;
        const double xhalf = 0.5 * x;

        y = y * (1.5 - fma(xhalf, y * y, 0.0));
        y = y * (1.5 - fma(xhalf, y * y, 0.0));

        return x * y;
    }

    // -------------------------------------------------------------------------
    // Phase-native inverse trig (returns turns)
    // -------------------------------------------------------------------------
    // These functions return angles in "turns" (phase):
    //   asin(z) ∈ [-0.25, +0.25],  acos(x) ∈ [0.0, 0.5],  atan(y) ∈ (-0.25, +0.25)
    // Internally they evaluate in radians (via atan2), then convert to phase.
    // They do not require or perform external phase wrapping helpers.

    static inline double asin(double z) {
        // ±π/2 in phase
        constexpr double HALFPI_PHASE = 0.25;  // (TAU * 0.25) / TAU
        if (z >= 1.0)  return HALFPI_PHASE;
        if (z <= -1.0) return -HALFPI_PHASE;

        // Near |z|≈1, use cancellation-safe approximation:
        // asin(z) ≈ sign(z) * (π/2 - sqrt(2*(1-|z|)))
        // Compute correction in radians, then convert to phase.
        double az = fabs(z);
        if (fabs(1.0 - az) < 1e-14) {
            double eps = 1.0 - az;
            double corr_rad = fast_sqrt(2.0 * eps);
            double corr_phase = corr_rad * INV_TAU;
            return (z >= 0.0 ? HALFPI_PHASE - corr_phase
                : -HALFPI_PHASE + corr_phase);
        }

        // asin(z) = atan2(z, sqrt(1 - z*z)) (radians) → convert to phase
        const double t = std::fma(-z, z, 1.0);          // 1 - z*z
        const double c = (t > 0.0) ? fast_sqrt(t) : 0.0;
        return atan2(z, c) * INV_TAU;
    }

    static inline double acos(double x) {
        // Clamp edges in phase units
        if (x >= 1.0)  return 0.0;       // 0 rad = 0 turns
        if (x <= -1.0) return 0.5;       // π rad = 1/2 turn

        // acos(x) = atan2( sqrt(1-x*x), x ) (radians) → convert to phase
        const double t = std::fma(-x, x, 1.0);          // 1 - x*x
        const double s = (t > 0.0) ? fast_sqrt(t) : 0.0;
        return atan2(s, x) * INV_TAU;
    }

    static __forceinline double atan(double y) {
        // atan(y) = atan2(y,1) (radians) → convert to phase
        return atan2(y, 1.0) * INV_TAU;
    }

} // namespace btrig

#endif // BTRIG_H
