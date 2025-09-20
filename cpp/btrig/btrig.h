#pragma once
#ifndef BTRIG_H
#define BTRIG_H

#include <cstdint>
#include <cfloat>
#include <cmath>

#if defined(_MSC_VER) && !defined(__clang__)
  #include <intrin.h>
#endif

namespace btrig {

    // =========================================================================
    // Q64.64 phase: canonical angle representation in turns [0,1)
    // NOTE:
    //  - On GCC/Clang we use 128-bit for fractions; ONE_Q64 is a true integer.
    //  - On MSVC (no __int128), fractions are 64-bit; ONE_Q64_L is used for scaling.
    // =========================================================================
#if defined(_MSC_VER) && !defined(__clang__)
    using q64  = long long;          // signed 64
    using uq64 = unsigned long long; // unsigned 64
    // 2^64 as long double for scaling (MSVC long double == double; still OK)
    constexpr long double ONE_Q64_L = 18446744073709551616.0L; // 2^64
#else
    using q64  = __int128_t;
    using uq64 = __uint128_t;
    constexpr uq64 ONE_Q64 = (uq64)1 << 64;
#endif

    constexpr int  FRAC_BITS = 64;

#if defined(_MSC_VER) && !defined(__clang__)
    static_assert(LDBL_MANT_DIG >= 53, "MSVC long double has 53-bit mantissa (same as double).");
#else
    static_assert(LDBL_MANT_DIG >= 64, "long double must have >=64-bit mantissa");
#endif

    // =========================================================================
    // Configuration: anchor density over one full turn (power of two)
    // =========================================================================
    constexpr int TRIG_ANCHORS_FULL_POW = 8;  // 256 anchors
    static_assert(TRIG_ANCHORS_FULL_POW >= 1, "bad TRIG_ANCHORS_FULL_POW");

    constexpr int NUM_ANCHORS_FULL = 1 << TRIG_ANCHORS_FULL_POW;
    constexpr int ANCHOR_MASK      = NUM_ANCHORS_FULL - 1;

    // =========================================================================
    // Core constants (long double)
    // =========================================================================
    constexpr long double TAU  = 6.283185307179586476925286766559005768394338798750211641949L;
    constexpr long double QUAD = TAU / 4.0L;                 // pi/2
    constexpr long double STEP = TAU / (long double)NUM_ANCHORS_FULL;

    // Reciprocal factorials (for tiny Horner polys)
    constexpr long double INV_2         = 1.0L / 2.0L;
    constexpr long double INV_6         = 1.0L / 6.0L;
    constexpr long double INV_24        = 1.0L / 24.0L;
    constexpr long double INV_120       = 1.0L / 120.0L;
    constexpr long double INV_720       = 1.0L / 720.0L;
    constexpr long double INV_5040      = 1.0L / 5040.0L;
    constexpr long double INV_40320     = 1.0L / 40320.0L;
    constexpr long double INV_362880    = 1.0L / 362880.0L;     // 9!
    constexpr long double INV_3628800   = 1.0L / 3628800.0L;    // 10!
    constexpr long double INV_39916800  = 1.0L / 39916800.0L;   // 11!
    constexpr long double INV_479001600 = 1.0L / 479001600.0L;  // 12!

    // =========================================================================
    // Anchor table: [i][0]=phase(turns), [i][1]=cos, [i][2]=sin
    // =========================================================================
    static long double fullCircleAnchors[NUM_ANCHORS_FULL][3];

    // =========================================================================
    // Portable helpers
    // =========================================================================
    // modfl shim (some libcs expose only ::modfl; MSVC prefers std::modf)
    static inline long double modf_ld(long double x, long double* ip) {
    #if defined(_MSC_VER) && !defined(__clang__)
        return std::modf(x, ip);
    #else
        return ::modfl(x, ip);
    #endif
    }

    // phase -> anchor index + 64-bit fractional remainder
    static inline void phaseQ64_to_idx_rem(uq64 phase_q64, int& idx, uint64_t& rem64) {
        const unsigned long long N = (unsigned long long)NUM_ANCHORS_FULL;
    #if defined(_MSC_VER) && !defined(__clang__)
        unsigned long long hi;
        unsigned long long lo = _umul128((unsigned long long)phase_q64, N, &hi);
        idx   = (int)hi & ANCHOR_MASK;
        rem64 = (uint64_t)lo;
    #else
        __uint128_t t = (__uint128_t)phase_q64 * ( (__uint128_t)N );
        idx   = (int)(t >> 64) & ANCHOR_MASK;
        rem64 = (uint64_t)(t & (((__uint128_t)1 << 64) - 1));
    #endif
    }

    static inline long double rem64_to_ld_unit(uint64_t r) {
        return (long double)r * ::ldexpl(1.0L, -64);
    }

    // Convert radians -> Q64.64 turns, round-to-nearest
    static inline uq64 radians_to_q64_turns(long double rad) {
        long double ip;
        long double turns = modf_ld(rad / TAU, &ip);
        if (turns < 0.0L) turns += 1.0L;
    #if defined(_MSC_VER) && !defined(__clang__)
        long double scaled = turns * ONE_Q64_L;
        if (scaled <= 0.0L) return 0;
        if (scaled >= ONE_Q64_L) return (uq64)0;
        return (uq64)(scaled + 0.5L);
    #else
        long double scaled = turns * (long double)ONE_Q64;
        if (scaled <= 0.0L) return 0;
        if (scaled >= (long double)ONE_Q64) return (uq64)0;
        return (uq64)(scaled + 0.5L);
    #endif
    }

    // Convert turns in [0,1) (long double) -> Q64.64 turns, round-to-nearest
    static inline uq64 turns_to_q64_rn(long double turns) {
        long double ip;
        long double frac = modf_ld(turns, &ip);
        if (frac < 0.0L) frac += 1.0L;
    #if defined(_MSC_VER) && !defined(__clang__)
        long double scaled = frac * ONE_Q64_L;
        if (scaled <= 0.0L) return 0;
        if (scaled >= ONE_Q64_L) return (uq64)0;
        return (uq64)(scaled + 0.5L);
    #else
        long double scaled = frac * (long double)ONE_Q64;
        if (scaled <= 0.0L) return 0;
        if (scaled >= (long double)ONE_Q64) return (uq64)0;
        return (uq64)(scaled + 0.5L);
    #endif
    }

    // Convert Q64.64 phase -> long double turns in [0,1)
    static inline long double q64_to_turns(btrig::uq64 p) {
    #if defined(_MSC_VER) && !defined(__clang__)
        // MSVC path: fraction stored in 64 bits
        const long double TWO64 = 18446744073709551616.0L; // 2^64
        return (long double)p / TWO64;
    #else
        // GCC/Clang: fraction is the low 64 bits of __uint128_t
        const __uint128_t mask = (((__uint128_t)1) << 64) - 1;
        uint64_t frac = (uint64_t)(p & mask);
        return (long double)frac * ldexpl(1.0L, -64);
    #endif
    }

    // Convenience: Q64.64 -> radians
    static inline long double q64_to_radians(btrig::uq64 p) {
        return q64_to_turns(p) * btrig::TAU;
    }


    constexpr long double INV_DEG = 1.0L / 360.0L;

    // Robust to very large |deg|:
    static inline btrig::uq64 from_degrees_q64(long double deg) {
        long double d = ::fmodl(deg, 360.0L);
        if (d < 0.0L) d += 360.0L;           // map to [0, 360)
        return turns_to_q64_rn(d * INV_DEG);
    }

    static inline long double to_degrees_from_q64(btrig::uq64 phase_q64) {
        return q64_to_turns(phase_q64) * 360.0L;
    }


    // =========================================================================
    // Helpers for anchors & polynomials
    // =========================================================================
    static inline void applyQuadrantTransform(long double inx, long double iny, int quadrant,
                                              long double& outx, long double& outy) {
        switch (quadrant & 3) {
        case 0: outx =  inx; outy =  iny; break;
        case 1: outx = -iny; outy =  inx; break;
        case 2: outx = -inx; outy = -iny; break;
        default:outx =  iny; outy = -inx; break;
        }
    }

    // Horner polynomials on [0, pi/2] for base anchors
    static inline long double cos_horner_ld(long double x, int max_even) {
        const long double z = x * x;
        int N = max_even / 2;
        long double acc = 1.0L;
        for (int k = N; k > 0; --k) {
            acc = 1.0L - (z * acc) / ((2*k) * (2*k - 1));
        }
        return acc;
    }
    static inline long double sin_horner_ld(long double x, int max_odd) {
        const long double z = x * x;
        int N = (max_odd - 1) / 2;
        long double acc = 1.0L;
        for (int k = N; k > 0; --k) {
            acc = 1.0L - (z * acc) / ((2*k + 1) * (2*k));
        }
        return x * acc;
    }

    // =========================================================================
    // Anchor initialization
    // =========================================================================
    static inline void initFullCircleAnchors() {
        constexpr int COS_MAX_EVEN = 22;
        constexpr int SIN_MAX_ODD  = 21;

        for (int i = 0; i < NUM_ANCHORS_FULL; ++i) {
            const long double phase  = (long double)i / (long double)NUM_ANCHORS_FULL;
            const long double angle  = phase * TAU;
            const int         quad   = (int)(angle / QUAD);
            const long double x      = angle - (long double)quad * QUAD;

            const long double c0 = cos_horner_ld(x, COS_MAX_EVEN);
            const long double s0 = sin_horner_ld(x, SIN_MAX_ODD);

            long double cA, sA;
            applyQuadrantTransform(c0, s0, quad, cA, sA);

            fullCircleAnchors[i][0] = phase;
            fullCircleAnchors[i][1] = cA;
            fullCircleAnchors[i][2] = sA;
        }

        // pin cardinals
        constexpr int Q = NUM_ANCHORS_FULL / 4;
        fullCircleAnchors[0][1]      =  1.0L; fullCircleAnchors[0][2]      =  0.0L;
        fullCircleAnchors[1 * Q][1]  =  0.0L; fullCircleAnchors[1 * Q][2]  =  1.0L;
        fullCircleAnchors[2 * Q][1]  = -1.0L; fullCircleAnchors[2 * Q][2]  =  0.0L;
        fullCircleAnchors[3 * Q][1]  =  0.0L; fullCircleAnchors[3 * Q][2]  = -1.0L;
    }

    static struct _InitFullCircleAnchors {
        _InitFullCircleAnchors() { initFullCircleAnchors(); }
    } _initFullCircleAnchorsInstance;

    // Residual scale: TAU / (N * 2^64)
    constexpr long double TWO64_L     =
        18446744073709551616.0L; // 2^64
    constexpr long double N_ANCHORS_L = (long double)NUM_ANCHORS_FULL;
    constexpr long double DENOM_RES_L = N_ANCHORS_L * TWO64_L;
    constexpr long double K_RESIDUAL  = TAU / DENOM_RES_L;

    // =========================================================================
    // Core: sincos from Q64.64 phase
    //  - precise=true  : cos(d) up to d^8/8!, sin(d) up to d^7/7!
    //  - precise=false : cos(d) ~ 1 - d^2/2 + d^4/24, sin(d) ~ d - d^3/6
    // =========================================================================
    static inline void sincos_q64(uq64 phase_q64, bool precise, long double& cosOut, long double& sinOut) {
        int idx; uint64_t rem;
        phaseQ64_to_idx_rem(phase_q64, idx, rem);

        // nearest-anchor pivot (right-edge ownership)
        const uint64_t sign = rem >> 63;
        const int idx2 = (idx + (int)sign) & ANCHOR_MASK;

        // centered remainder -> radians
        const int64_t srem = (int64_t)rem;
        const long double d = (long double)srem * K_RESIDUAL;
        const long double z = d * d;

        const long double* a = fullCircleAnchors[idx2];
        const long double cA = a[1], sA = a[2];

        long double dx, dy;
        long double ts, tc;
        if (precise) {
            // cos(d) = 1 - z/2! + z^2/4! - z^3/6! + z^4/8!
            tc =  -INV_40320;   // -1/8!
            tc = tc * z +  INV_720;         // +1/6!
            tc = tc * z -  INV_24;          // -1/4!
            tc = tc * z +  INV_2;           // +1/2!
            dx = 1.0L - z * tc;

            // sin(d) = d * (1 - z/3! + z^2/5! - z^3/7!)
            ts =  -INV_5040;    // -1/7!
            ts = ts * z +  INV_120;         // +1/5!
            ts = ts * z -  INV_6;           // -1/3!
            ts = ts * z +  1.0L;            // +1
            dy = d * ts;
        } else {
            // fast 4/3-term
            tc = -INV_24;       // -1/4!
            tc = tc * z +  INV_2;           // +1/2!
            dx = 1.0L - z * tc;

            ts = -INV_6;        // -1/3!
            ts = ts * z +  1.0L;            // +1
            dy = d * ts;
        }

        // rotate anchor by residual
        cosOut = dx * cA - dy * sA;
        sinOut = dx * sA + dy * cA;
    }

    // --- internal micro-rotation (residual) kernel --------------------
    static inline void _btrig_residual(long double d, bool precise,
                                    long double& dx, long double& dy) {
        const long double z = d * d;
        long double tc, ts;

        if (precise) {
            // cos(d) = 1 - z/2! + z^2/4! - z^3/6! + z^4/8!
            tc =  -btrig::INV_40320;    // -1/8!
            tc = tc * z +  btrig::INV_720;          // +1/6!
            tc = tc * z -  btrig::INV_24;           // -1/4!
            tc = tc * z +  btrig::INV_2;            // +1/2!
            dx = 1.0L - z * tc;

            // sin(d) = d * (1 - z/3! + z^2/5! - z^3/7!)
            ts =  -btrig::INV_5040;     // -1/7!
            ts = ts * z +  btrig::INV_120;          // +1/5!
            ts = ts * z -  btrig::INV_6;            // -1/3!
            ts = ts * z +  1.0L;                    // +1
            dy = d * ts;
        } else {
            // fast: cos(d) ~ 1 - d^2/2 + d^4/24 ; sin(d) ~ d - d^3/6
            tc = -btrig::INV_24;        // -1/4!
            tc = tc * z +  btrig::INV_2;            // +1/2!
            dx = 1.0L - z * tc;

            ts = -btrig::INV_6;         // -1/3!
            ts = ts * z +  1.0L;                    // +1
            dy = d * ts;
        }
    }

    // --- cosine(single) ----------------------------------------------
    static inline long double cos_q64(btrig::uq64 phase_q64, bool precise) {
        int idx; uint64_t rem;
        btrig::phaseQ64_to_idx_rem(phase_q64, idx, rem);

        // nearest-anchor pivot (right-edge ownership)
        const uint64_t sign = rem >> 63;
        const int idx2 = (idx + (int)sign) & btrig::ANCHOR_MASK;

        // centered remainder -> radians
        const int64_t srem = (int64_t)rem;
        const long double d = (long double)srem * btrig::K_RESIDUAL;

        const long double* a = btrig::fullCircleAnchors[idx2];
        const long double cA = a[1], sA = a[2];

        long double dx, dy;
        _btrig_residual(d, precise, dx, dy);

        // cos(x+d) = dx*cA - dy*sA
        return dx * cA - dy * sA;
    }

    // --- sine(single) -------------------------------------------------
    static inline long double sin_q64(btrig::uq64 phase_q64, bool precise) {
        int idx; uint64_t rem;
        btrig::phaseQ64_to_idx_rem(phase_q64, idx, rem);

        // nearest-anchor pivot (right-edge ownership)
        const uint64_t sign = rem >> 63;
        const int idx2 = (idx + (int)sign) & btrig::ANCHOR_MASK;

        // centered remainder -> radians
        const int64_t srem = (int64_t)rem;
        const long double d = (long double)srem * btrig::K_RESIDUAL;

        const long double* a = btrig::fullCircleAnchors[idx2];
        const long double cA = a[1], sA = a[2];

        long double dx, dy;
        _btrig_residual(d, precise, dx, dy);

        // sin(x+d) = dx*sA + dy*cA
        return dx * sA + dy * cA;
    }

    // =========================================================================
    // Tangent from Q64.64 phase (computed as sin/cos of rotated anchor)
    // =========================================================================
    static inline long double tan_q64(uq64 phase_q64, bool precise) {
        long double c, s;
        sincos_q64(phase_q64, precise, c, s);
        return s / c; // naturally tends to +/-inf near verticals
    }

    // =========================================================================
    // atan / atan2 core (radians), robust argument reduction
    // =========================================================================
    // atan polynomial on small |u| (|u| <= ~0.4142) : choose degree ∈ {3,5,7,9}
    static inline long double atan_poly_small(long double u, int deg) {
        const long double z = u*u;
        // x - x^3/3 + x^5/5 - x^7/7 + x^9/9
        constexpr long double c1 =  1.0L;
        constexpr long double c3 = -1.0L/3.0L;
        constexpr long double c5 =  1.0L/5.0L;
        constexpr long double c7 = -1.0L/7.0L;
        constexpr long double c9 =  1.0L/9.0L;
        long double acc;
        switch (deg) {
            case 9:  acc = c9; acc = c7 + z*acc; acc = c5 + z*acc; acc = c3 + z*acc; break;
            case 7:  acc = c7; acc = c5 + z*acc; acc = c3 + z*acc;                    break;
            case 5:  acc = c5; acc = c3 + z*acc;                                      break;
            default: acc = c3;                                                        break; // deg=3
        }
        return u * (c1 + z*acc);
    }

    // Reduce x >= 0 using:
    //  - if x > 1      : atan(x) = pi/2 - atan(1/x)
    //  - if x > tan(pi/8): atan(x) = pi/4 + atan((x-1)/(1+x))
    static inline long double atan_reduce(long double x, bool precise) {
        const long double T = 0.4142135623730950488016887242096980785696718753769480732L; // tan(pi/8)
        int q = 0; // quarter-step flags: bit1=pi/2, bit0=pi/4
        long double u = x;

        if (u > 1.0L) { u = 1.0L/u; q ^= 2; }                  // add pi/2 later
        if (u > T)    { u = (u - 1.0L) / (1.0L + u); q ^= 1; } // add pi/4 later

        long double a = precise ? atan_poly_small(u, 9)
                                : (u + (u*u*u)/3.0L);          // x + x^3/3

        if (q & 1) a += QUAD / 2.0L; // +pi/4
        if (q & 2) a  = QUAD - a;    // pi/2 - a
        return a;
    }

    // atan2 in radians, return in [-pi, pi]
    static inline long double atan2_ld(long double y, long double x, bool precise) {
        if (x > 0.0L) {
            long double a = atan_reduce(::fabsl(y/x), precise);
            return ::copysignl(a, y);
        } else if (x < 0.0L) {
            long double a = atan_reduce(::fabsl(y/x), precise);
            long double base = TAU / 2.0L; // pi
            return (y >= 0.0L) ? (base - a) : (-base + a);
        } else { // x == 0
            return (y > 0.0L) ? QUAD : (y < 0.0L ? -QUAD : 0.0L);
        }
    }

    // =========================================================================
    // Public angle → Q64.64 converters
    //  - atan_q64: angle in [-pi/2, pi/2] mapped to turns in [−0.25, 0.25) mod 1
    //  - atan2_q64: full-circle angle
    //  - asin_q64, acos_q64 via stable atan2 forms
    // =========================================================================
    static inline uq64 atan_q64(long double x, bool precise) {
        long double ang = atan_reduce(::fabsl(x), precise);
        ang = ::copysignl(ang, x);
        return radians_to_q64_turns(ang);
    }

    static inline uq64 atan2_q64(long double y, long double x, bool precise) {
        return radians_to_q64_turns(atan2_ld(y, x, precise));
    }

    static inline uq64 asin_q64(long double x, bool precise) {
        if (x >  1.0L) x =  1.0L;
        if (x < -1.0L) x = -1.0L;
        long double y  = x;
        long double xr = ::sqrtl(::fmaxl(0.0L, 1.0L - x*x));
        long double ang = atan2_ld(y, xr, precise);   // [-pi/2, pi/2]
        return radians_to_q64_turns(ang);
    }

    static inline uq64 acos_q64(long double x, bool precise) {
        if (x >  1.0L) x =  1.0L;
        if (x < -1.0L) x = -1.0L;
        long double y  = ::sqrtl(::fmaxl(0.0L, 1.0L - x*x));
        long double ang = atan2_ld(y, x, precise);    // [0, pi]
        return radians_to_q64_turns(ang);
    }

    #if defined(_MSC_VER) && !defined(__clang__)
    static inline long double fast_sqrt_ld_msvc(long double x) {
        if (!(x > 0.0L)) return (x == 0.0L) ? 0.0L : std::numeric_limits<long double>::quiet_NaN();
        // Direct double version: two NR steps for good accuracy
        double xd = (double)x;
        uint64_t bits; std::memcpy(&bits, &xd, 8);
        bits = 0x5fe6ec85e7de30daULL - (bits >> 1);
        double y; std::memcpy(&y, &bits, 8);
        y = y * (1.5 - 0.5 * xd * y * y);
        y = y * (1.5 - 0.5 * xd * y * y);
        return (long double)(xd * y);
    }
#endif


    // Constants
    constexpr long double LN2    = 0.6931471805599453094172321214581765680755L;
    constexpr long double INV_LN2= 1.4426950408889634073599246810018921374266L;

    // exp core on small |r| (r in approx [-LN2/2, LN2/2])
    // Horner form; precise=true uses degree 10, false uses degree 6.
    static inline long double _exp_poly(long double r, bool precise) {
        const long double r2 = r * r;
        if (precise) {
            // 1 + r + r^2/2! + ... + r^10/10!
            // Horner from highest degree for stability.
            long double p = 1.0L/3628800.0L;     // 1/10!
            p = p * r + 1.0L/362880.0L;          // 1/9!
            p = p * r + 1.0L/40320.0L;           // 1/8!
            p = p * r + 1.0L/5040.0L;            // 1/7!
            p = p * r + 1.0L/720.0L;             // 1/6!
            p = p * r + 1.0L/120.0L;             // 1/5!
            p = p * r + 1.0L/24.0L;              // 1/4!
            p = p * r + 1.0L/6.0L;               // 1/3!
            p = p * r + 1.0L/2.0L;               // 1/2!
            p = p * r + 1.0L;                    // 1
            return 1.0L + r * p;
        } else {
            // degree-6: up through r^6/6!
            long double p = 1.0L/720.0L;         // 1/6!
            p = p * r + 1.0L/120.0L;             // 1/5!
            p = p * r + 1.0L/24.0L;              // 1/4!
            p = p * r + 1.0L/6.0L;               // 1/3!
            p = p * r + 1.0L/2.0L;               // 1/2!
            p = p * r + 1.0L;                    // 1
            return 1.0L + r * p;
        }
    }

    // exp for long double
    static inline long double exp_ld(long double x, bool precise) {
        if (x != x) return x;                // NaN
        if (x > 11356.0L) return INFINITY;   // huge -> inf (safe guard)
        if (x < -11356.0L) return 0.0L;      // tiny -> 0

        // Range reduction: x = k*LN2 + r, with r small
        long double kf = nearbyintl(x * INV_LN2); // k as long double integer value
        long double r  = x - kf * LN2;

        // Core
        long double er = _exp_poly(r, precise);

        // Scale back by 2^k
        // k might be outside int range on some platforms; convert carefully:
        long long k = (long long)kf;
        return ldexpl(er, (int)k);
    }

    // log core: log(m) with m in (0, +inf)
    // We map m to y=(m-1)/(m+1) so |y| <= ~0.1716 over m in [sqrt(0.5), sqrt(2))
    // log(m) = 2 * atanh(y) ≈ 2 * (y + y^3/3 + y^5/5 + ...)
    static inline long double _log_core_atanh(long double m, bool precise) {
        // y in (-0.1716..0.1716)
        long double y  = (m - 1.0L) / (m + 1.0L);
        long double z  = y * y;

        if (precise) {
            // up to y^11/11
            long double p = 1.0L/11.0L;
            p = p * z + 1.0L/9.0L;
            p = p * z + 1.0L/7.0L;
            p = p * z + 1.0L/5.0L;
            p = p * z + 1.0L/3.0L;
            return 2.0L * (y + y * z * p);
        } else {
            // up to y^7/7
            long double p = 1.0L/7.0L;
            p = p * z + 1.0L/5.0L;
            p = p * z + 1.0L/3.0L;
            return 2.0L * (y + y * z * p);
        }
    }

    // log for long double (natural log)
    static inline long double log_ld(long double x, bool precise) {
        if (x < 0.0L)  return NAN;
        if (x == 0.0L) return -INFINITY;
        if (!std::isfinite(x)) return x; // +inf or NaN

        // Decompose: x = m * 2^e with m in [0.5,1)
        int e;
        long double m = frexpl(x, &e);

        // Normalize m into [sqrt(0.5), sqrt(2)) for best core accuracy
        // sqrt(0.5) ~= 0.70710678, sqrt(2) ~= 1.41421356
        if (m < 0.70710678118654752440L) {
            m *= 2.0L;
            e -= 1;
        }

        long double core = _log_core_atanh(m, precise);
        // Combine with exponent term: e * ln(2)
        return core + (long double)e * LN2;
    }

    // Optional double convenience wrappers
    static inline double exp_fast(double x)   { return (double)exp_ld((long double)x, false); }
    static inline double exp_precise(double x){ return (double)exp_ld((long double)x, true ); }
    static inline double log_fast(double x)   { return (double)log_ld((long double)x, false); }
    static inline double log_precise(double x){ return (double)log_ld((long double)x, true ); }


} // namespace btrig

#endif // BTRIG_H
