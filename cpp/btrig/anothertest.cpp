// btrig_tests_q64.cpp
#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>
#include <algorithm>
#include "btrig.h"

// ---------- portable helpers ----------
#if defined(_MSC_VER) && !defined(__clang__)
    #define BTRIG_MSVC 1
#else
    #define BTRIG_MSVC 0
#endif

struct Stats {
    long double sumsq = 0.0L;
    long double maxabs = 0.0L;
    uint64_t    count = 0;
    void add(long double e) {
        long double a = fabsl(e);
        sumsq += a * a;
        if (a > maxabs) maxabs = a;
        ++count;
    }
    long double rms() const { return (count ? sqrtl(sumsq / (long double)count) : 0.0L); }
};

static inline long double q64_to_turns(btrig::uq64 p) {
#if BTRIG_MSVC
    // MSVC stores fractional turns in 64-bit; scale by 2^-64
    const long double TWO64 = 18446744073709551616.0L; // 2^64
    return (long double)p / TWO64;
#else
    // GCC/Clang use __uint128_t; fraction is in the low 64 bits
    const __uint128_t mask = (((__uint128_t)1) << 64) - 1;
    uint64_t low = (uint64_t)(p & mask);
    return (long double)low * ldexpl(1.0L, -64);
#endif
}
static inline long double q64_to_radians(btrig::uq64 p) {
    return q64_to_turns(p) * btrig::TAU;
}
static inline long double wrap_pm_pi(long double a) {
    // map to (-pi, pi]
    long double t = fmodl(a, btrig::TAU);
    if (t <= -btrig::TAU * 0.5L) t += btrig::TAU;
    if (t >   btrig::TAU * 0.5L) t -= btrig::TAU;
    return t;
}
static inline long double angdiff(long double a, long double b) {
    return fabsl(wrap_pm_pi(a - b));
}

static inline btrig::uq64 random_phase_q64(std::mt19937_64& rng) {
#if BTRIG_MSVC
    return (btrig::uq64)rng();
#else
    __uint128_t hi = (__uint128_t)rng();
    __uint128_t lo = (__uint128_t)rng();
    return (hi << 64) | lo;
#endif
}

// ---------- tests ----------
void test_anchors_vs_libm() {
    using namespace btrig;
    Stats sc, ss;
    for (int i = 0; i < NUM_ANCHORS_FULL; ++i) {
        long double phase = (long double)i / (long double)NUM_ANCHORS_FULL;
        long double ang   = phase * TAU;
        long double c_std = cosl(ang);
        long double s_std = sinl(ang);

        long double c_a = fullCircleAnchors[i][1];
        long double s_a = fullCircleAnchors[i][2];

        sc.add(c_a - c_std);
        ss.add(s_a - s_std);
    }
    std::cout << "=== Anchors vs libm ===\n";
    std::cout << std::setprecision(18)
              << "rms|dcos(anchor)|=" << (double)sc.rms()
              << "   max|dcos(anchor)|=" << (double)sc.maxabs << "\n"
              << "rms|dsin(anchor)|=" << (double)ss.rms()
              << "   max|dsin(anchor)|=" << (double)ss.maxabs << "\n\n";
}

void test_random_sincos(bool precise, int N, uint64_t seed = 0xBEEFFACE1234ULL) {
    using namespace btrig;
    std::mt19937_64 rng(seed);
    Stats ec, es;

    for (int i = 0; i < N; ++i) {
        uq64 p = random_phase_q64(rng);

        long double c, s;
        sincos_q64(p, precise, c, s);

        long double ang = q64_to_radians(p);
        long double c_std = cosl(ang);
        long double s_std = sinl(ang);

        ec.add(c - c_std);
        es.add(s - s_std);
    }

    std::cout << "=== Random phases, precise=" << (precise ? "true" : "false")
              << " (N=" << N << ") ===\n";
    std::cout << std::setprecision(18)
              << "rms|dcos|=" << (double)ec.rms()
              << "   max|dcos|=" << (double)ec.maxabs << "\n"
              << "rms|dsin|=" << (double)es.rms()
              << "   max|dsin|=" << (double)es.maxabs << "\n\n";
}

void test_boundary_probes(bool precise) {
    using namespace btrig;
    std::vector<long double> centers = {
        0.00L, 0.25L, 0.50L, 0.75L, // cardinals
        0.125L, 0.375L, 0.625L, 0.875L // diagonals
    };
    // tiny steps around centers (in turns)
    std::vector<long double> deltas = {
        -1.0L/ (1ULL<<40), -1.0L/ (1ULL<<44), -1.0L/ (1ULL<<48),
         0.0L,
         1.0L/ (1ULL<<48),  1.0L/ (1ULL<<44),  1.0L/ (1ULL<<40)
    };

    Stats ec, es;
    for (long double center : centers) {
        for (long double d : deltas) {
            long double t = center + d;
            // normalize to [0,1)
            long double ip;
            t = modfl(t, &ip);
            if (t < 0.0L) t += 1.0L;

            btrig::uq64 p = btrig::turns_to_q64_rn(t);

            long double c, s;
            btrig::sincos_q64(p, precise, c, s);

            long double ang = q64_to_radians(p);
            long double c_std = cosl(ang);
            long double s_std = sinl(ang);

            ec.add(c - c_std);
            es.add(s - s_std);
        }
    }
    std::cout << "=== Boundary probes (cardinals + diagonals), precise=" << (precise?"true":"false") << " ===\n";
    std::cout << std::setprecision(18)
              << "rms|dcos|=" << (double)ec.rms()
              << "   max|dcos|=" << (double)ec.maxabs << "\n"
              << "rms|dsin|=" << (double)es.rms()
              << "   max|dsin|=" << (double)es.maxabs << "\n\n";
}

void test_tan(int N = 20000, uint64_t seed = 0xA11CE5EEDULL) {
    using namespace btrig;
    std::mt19937_64 rng(seed);
    Stats et;
    int skipped = 0;

    for (int i = 0; i < N; ++i) {
        uq64 p = random_phase_q64(rng);
        long double c, s;
        sincos_q64(p, true, c, s);

        // Skip near vertical to avoid inf/NaN comparisons
        if (fabsl(c) < 1e-30L) { skipped++; continue; }

        long double t_b = s / c;
        long double ang = q64_to_radians(p);
        long double t_std = tanl(ang);

        et.add(t_b - t_std);
    }
    std::cout << "=== tan_q64 via sin/cos (precise=true), random phases N=" << N
              << " (skipped " << skipped << " near-vertical) ===\n";
    std::cout << std::setprecision(18)
              << "rms|dtan|=" << (double)et.rms()
              << "   max|dtan|=" << (double)et.maxabs << "\n\n";
}

void test_inverse_trig(int N = 20000, uint64_t seed = 0xC0FFEEBADD00DULL) {
    using namespace btrig;
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<long double> uni(-10.0L, 10.0L);
    std::uniform_real_distribution<long double> unit(-1.0L, 1.0L);

    Stats ea, ea2, easin, eacos;

    // atan
    for (int i = 0; i < N; ++i) {
        long double x = uni(rng);
        uq64 t_q = atan_q64(x, /*precise=*/true);
        long double a_b = q64_to_radians(t_q);
        if (a_b >= btrig::TAU * 0.5L) a_b -= btrig::TAU; // map to [-pi, pi)

        long double a_std = atanl(x);
        ea.add(angdiff(a_b, a_std));
    }

    // atan2
    for (int i = 0; i < N; ++i) {
        long double y = uni(rng);
        long double x = uni(rng);
        uq64 t_q = atan2_q64(y, x, /*precise=*/true);
        long double a_b = q64_to_radians(t_q);
        if (a_b >= btrig::TAU * 0.5L) a_b -= btrig::TAU; // map to [-pi, pi)

        long double a_std = atan2l(y, x);
        ea2.add(angdiff(a_b, a_std));
    }

    // asin / acos
    for (int i = 0; i < N; ++i) {
        long double x = unit(rng);
        uq64 as_q = asin_q64(x, /*precise=*/true);
        long double as_b = q64_to_radians(as_q);
        if (as_b >= btrig::TAU * 0.5L) as_b -= btrig::TAU; // map to [-pi, pi)
        long double as_std = asinl(x);
        easin.add(angdiff(as_b, as_std));

        uq64 ac_q = acos_q64(x, /*precise=*/true);
        long double ac_b = q64_to_radians(ac_q);
        // acos is [0,pi], our mapping leaves it in [-pi,pi), which is fine
        long double ac_std = acosl(x);
        eacos.add(angdiff(ac_b, ac_std));
    }

    std::cout << "=== inverse trig (precise=true), random N=" << N << " ===\n";
    std::cout << std::setprecision(18)
              << "rms|datan|="  << (double)ea.rms()   << "   max|datan|="  << (double)ea.maxabs   << "\n"
              << "rms|datan2|=" << (double)ea2.rms()  << "   max|datan2|=" << (double)ea2.maxabs  << "\n"
              << "rms|dasin|="  << (double)easin.rms()<< "   max|dasin|="  << (double)easin.maxabs<< "\n"
              << "rms|dacos|="  << (double)eacos.rms()<< "   max|dacos|="  << (double)eacos.maxabs<< "\n\n";
}

template<class F>
static inline double time_ms_bestof3(F&& fn) {
    using clock = std::chrono::high_resolution_clock;
    double best = std::numeric_limits<double>::infinity();
    for (int k = 0; k < 3; ++k) {
        auto t0 = clock::now();
        fn();
        auto t1 = clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        if (ms < best) best = ms;
    }
    return best;
}

void bench_throughput(int N = 5000000) {
    using namespace btrig;
    std::mt19937_64 rng(0xDEADBEEFCAFEBABEULL);
    std::vector<uq64> phases(N);
    for (int i = 0; i < N; ++i) phases[i] = random_phase_q64(rng);

    volatile long double sinkC = 0.0L, sinkS = 0.0L; // prevent dead-code elimination

    double t_btrig = time_ms_bestof3([&](){
        long double c, s;
        for (int i = 0; i < N; ++i) {
            sincos_q64(phases[i], /*precise=*/false, c, s);
            sinkC += c; sinkS += s;
        }
    });

    double t_libm = time_ms_bestof3([&](){
        for (int i = 0; i < N; ++i) {
            long double ang = q64_to_radians(phases[i]);
            sinkC += cosl(ang);
            sinkS += sinl(ang);
        }
    });

    std::cout << "=== Throughput (best of 3), N=" << N << " ===\n"
              << std::fixed << std::setprecision(3)
              << "btrig::sincos_q64 (precise=false): " << t_btrig << " ms\n"
              << "libm sinl/cosl:                    " << t_libm  << " ms\n\n";
    (void)sinkC; (void)sinkS;
}

int main() {
    std::cout.setf(std::ios::scientific);
    std::cout << std::setprecision(18);

    test_anchors_vs_libm();
    test_random_sincos(/*precise=*/true,  20000);
    test_random_sincos(/*precise=*/false, 20000);
    test_boundary_probes(/*precise=*/true);
    test_tan();
    test_inverse_trig();
    bench_throughput(5000000);

    return 0;
}
