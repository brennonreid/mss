// test_btrig_q64.cpp
#include <chrono>
#include <iostream>
#include <iomanip>
#include <random>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>
#include "btrig.h"

using namespace btrig;

// ------------------------------------------------------------------
// Compiler barriers to prevent dead-code elimination / constant fold
// ------------------------------------------------------------------
#if defined(__GNUC__) || defined(__clang__)
  #define NOINLINE __attribute__((noinline))
  #define HOT      __attribute__((hot))
#else
  #define NOINLINE
  #define HOT
#endif

static inline void ClobberMemory() {
#if defined(__GNUC__) || defined(__clang__)
  asm volatile("" : : : "memory");
#endif
}

template <class T>
static inline void DoNotOptimize(T const& v) {
#if defined(__GNUC__) || defined(__clang__)
  asm volatile("" : : "g"(v) : "memory");
#endif
}

// Global volatile sink so results must be materialized.
static volatile long double g_sink = 0.0L;

// ------------------------------------------------------------------
// Portable sincosl / sincos shims + volatile function pointers
// ------------------------------------------------------------------
static inline void sincosl_portable(long double a, long double* s, long double* c) {
#if defined(__GLIBC__) || defined(__APPLE__) || defined(__linux__)
    ::sincosl(a, s, c);
#else
    *s = ::sinl(a);
    *c = ::cosl(a);
#endif
}
static inline void sincos_portable(double a, double* s, double* c) {
#if defined(__GLIBC__) || defined(__APPLE__) || defined(__linux__)
    ::sincos(a, s, c);
#else
    *s = ::sin(a);
    *c = ::cos(a);
#endif
}
using sincosl_fn = void(*)(long double, long double*, long double*);
using sincos_fn  = void(*)(double,       double*,       double*);
static volatile sincosl_fn p_sincosl = &sincosl_portable;
static volatile sincos_fn  p_sincos  = &sincos_portable;

// ------------------------------------------------------------------
// Helpers
// ------------------------------------------------------------------

// Random Q64.64 in [0,1)
static inline std::vector<btrig::uq64> make_random_q64(size_t N, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<uint64_t> U64(0ULL, ~0ULL);
    std::vector<btrig::uq64> v(N);
    for (size_t i=0;i<N;++i) v[i] = (btrig::uq64)U64(rng);
    return v;
}

// Convert Q64.64 turns -> long double radians (for libm evals only)
static inline long double q64_to_radians(btrig::uq64 q) {
    long double frac = (long double)((long double)q * ldexpl(1.0L, -64)); // exact -> [0,1)
    return frac * btrig::TAU;
}

// Build btrig precise reference arrays
static inline void compute_btrig_ref(const std::vector<btrig::uq64>& phases_q,
                                     std::vector<long double>& refC,
                                     std::vector<long double>& refS)
{
    const size_t N = phases_q.size();
    refC.resize(N);
    refS.resize(N);
    for (size_t i=0;i<N;++i) {
        long double c,s;
        btrig::sincos_q64(phases_q[i], true, c, s);
        refC[i] = c;
        refS[i] = s;
    }
}

// Phase index -> exact Q64.64 (i/N, rounded)
static inline btrig::uq64 phase_index_to_q64(int i) {
    __uint128_t num = (__uint128_t)( (unsigned long long)i ) << 64;
    unsigned long long N = (unsigned long long)NUM_ANCHORS_FULL;
    __uint128_t q = num / N;
    __uint128_t r = num % N;
    if (r * 2 >= N) ++q;
    return (btrig::uq64)q;
}

/*
// Turns (long double) -> Q64.64 rounding
static inline btrig::uq64 turns_to_q64_rn(long double turns) {
    long double ip;
    long double frac = ::modfl(turns, &ip);
    if (frac < 0.0L) frac += 1.0L;
    long double scaled = frac * (long double)btrig::ONE_Q64;
    if (scaled <= 0.0L) return 0;
    if (scaled >= (long double)btrig::ONE_Q64) return (btrig::uq64)0;
    return (btrig::uq64)(scaled + 0.5L);
}
*/

/*
constexpr long double INV_DEG = 1.0L / 360.0L;
// from_degrees: long double -> uq64
static inline btrig::uq64 from_degrees_q64(long double deg) {
    // 360 degrees = 1 turn
    long double turns = deg * INV_DEG;
    return turns_to_q64_rn(turns);
}
*/


// Error accumulator vs reference arrays
struct ErrAgg {
    long double rmsC=0, rmsS=0, maxC=0, maxS=0;
    void add(long double c, long double s, long double rc, long double rs) {
        long double ec = fabsl(c - rc);
        long double es = fabsl(s - rs);
        rmsC += ec*ec; rmsS += es*es;
        if (ec > maxC) maxC = ec;
        if (es > maxS) maxS = es;
    }
    void finalize(size_t N) { rmsC = std::sqrt(rmsC/N); rmsS = std::sqrt(rmsS/N); }
};

// ------------------------------------------------------------------
// Results struct
// ------------------------------------------------------------------
struct BenchResult {
    const char* name;
    double secs;
    double ns_per_call;
    double mevals_per_s;
    long double rms_err_cos, rms_err_sin;
    long double max_err_cos, max_err_sin;
    long double sink;
};

// ------------------------------------------------------------------
// Bench: btrig::sincos_q64 (precise or fast), errors vs btrig precise ref
// ------------------------------------------------------------------
template<bool PRECISE>
NOINLINE static BenchResult bench_btrig_q64(const std::vector<btrig::uq64>& phases_q,
                                            const std::vector<long double>& refC,
                                            const std::vector<long double>& refS)
{
    using namespace std::chrono;
    const int N = (int)phases_q.size();

    // Warm-up
    {
        long double c,s;
        for (int i=0;i<1000;++i) btrig::sincos_q64(phases_q[i % N], PRECISE, c, s);
    }

    // Timed kernel
    long double sink = 0.0L;
    auto t0 = high_resolution_clock::now();
    for (int i=0;i<N;++i) {
        long double c,s;
        btrig::sincos_q64(phases_q[i], PRECISE, c, s);
        sink += c + s;
        DoNotOptimize(c);
        DoNotOptimize(s);
        ClobberMemory();
    }
    DoNotOptimize(sink);
    ClobberMemory();
    auto t1 = high_resolution_clock::now();

    // Accuracy vs btrig precise reference (not timed)
    ErrAgg agg;
    for (int i=0;i<N;++i) {
        long double c,s;
        btrig::sincos_q64(phases_q[i], PRECISE, c, s);
        agg.add(c, s, refC[i], refS[i]);
    }
    agg.finalize(N);

    g_sink += sink;
    ClobberMemory();

    double secs = duration<double>(t1 - t0).count();
    return { PRECISE ? "btrig::sincos_q64 (precise=true, ref=btrig precise)"
                     : "btrig::sincos_q64 (precise=false, vs btrig precise)",
             secs, secs*1e9/N, N/(secs*1e6), agg.rmsC, agg.rmsS, agg.maxC, agg.maxS, sink };
}

// ------------------------------------------------------------------
// Bench: libm sincos (double), errors vs btrig precise ref
// ------------------------------------------------------------------
NOINLINE static BenchResult bench_sincos_vs_btrig(const std::vector<btrig::uq64>& phases_q,
                                                  const std::vector<long double>& refC,
                                                  const std::vector<long double>& refS)
{
    using namespace std::chrono;
    const int N = (int)phases_q.size();

    std::vector<double> ang(N);
    for (int i=0;i<N;++i) ang[i] = (double)q64_to_radians(phases_q[i]);

    double s,c;
    for (int i=0;i<1000;++i) p_sincos(ang[i%N], &s, &c);

    long double sink = 0.0L;
    auto t0 = high_resolution_clock::now();
    for (int i=0;i<N;++i) {
        p_sincos(ang[i], &s, &c);
        sink += (long double)s + (long double)c;
        DoNotOptimize(s);
        DoNotOptimize(c);
        ClobberMemory();
    }
    DoNotOptimize(sink);
    ClobberMemory();
    auto t1 = high_resolution_clock::now();

    ErrAgg agg;
    for (int i=0;i<N;++i) {
        double sd, cd;
        p_sincos(ang[i], &sd, &cd);
        agg.add((long double)cd, (long double)sd, refC[i], refS[i]);
    }
    agg.finalize(N);

    g_sink += sink;
    ClobberMemory();

    double secs = duration<double>(t1 - t0).count();
    return {"libm sincos (double, vs btrig precise)", secs, secs*1e9/N, N/(secs*1e6),
            agg.rmsC, agg.rmsS, agg.maxC, agg.maxS, sink};
}

// ------------------------------------------------------------------
// Bench: libm sincosl (long double), errors vs btrig precise ref
// ------------------------------------------------------------------
NOINLINE static BenchResult bench_sincosl_vs_btrig(const std::vector<btrig::uq64>& phases_q,
                                                   const std::vector<long double>& refC,
                                                   const std::vector<long double>& refS)
{
    using namespace std::chrono;
    const int N = (int)phases_q.size();

    std::vector<long double> ang(N);
    for (int i=0;i<N;++i) ang[i] = q64_to_radians(phases_q[i]);

    long double s,c;
    for (int i=0;i<1000;++i) p_sincosl(ang[i%N], &s, &c);

    long double sink = 0.0L;
    auto t0 = high_resolution_clock::now();
    for (int i=0;i<N;++i) {
        p_sincosl(ang[i], &s, &c);
        sink += s + c;
        DoNotOptimize(s);
        DoNotOptimize(c);
        ClobberMemory();
    }
    DoNotOptimize(sink);
    ClobberMemory();
    auto t1 = high_resolution_clock::now();

    ErrAgg agg;
    for (int i=0;i<N;++i) {
        long double sd, cd;
        p_sincosl(ang[i], &sd, &cd);
        agg.add((long double)cd, (long double)sd, refC[i], refS[i]);
    }
    agg.finalize(N);

    g_sink += sink;
    ClobberMemory();

    double secs = duration<double>(t1 - t0).count();
    return {"libm sincosl (vs btrig precise)", secs, secs*1e9/N, N/(secs*1e6),
            agg.rmsC, agg.rmsS, agg.maxC, agg.maxS, sink};
}

NOINLINE static BenchResult bench_tan_q64(const std::vector<btrig::uq64>& phases_q) {
    using namespace std::chrono;
    const int N = (int)phases_q.size();
    long double sink = 0.0L;
    auto t0 = high_resolution_clock::now();
    for (int i = 0; i < N; ++i) {
        long double t = btrig::tan_q64(phases_q[i], true);
        sink += t;
        DoNotOptimize(t);
        ClobberMemory();
    }
    auto t1 = high_resolution_clock::now();
    g_sink += sink;
    double secs = duration<double>(t1 - t0).count();
    return {"btrig::tan_q64 (precise)", secs, secs*1e9/N, N/(secs*1e6), 0.0L, 0.0L, 0.0L, 0.0L, sink};
}

// ------------------------------------------------------------------
// Anchors: error vs btrig precise at anchors
// ------------------------------------------------------------------
static void test_anchors_vs_btrig() {
    std::cout << "=== Anchors: error vs btrig::sincos_q64(precise=true) ===\n";

    std::vector<btrig::uq64> phases;
    phases.reserve(NUM_ANCHORS_FULL);
    for (int i=0;i<NUM_ANCHORS_FULL;++i) phases.push_back(phase_index_to_q64(i));

    std::vector<long double> refC, refS;
    compute_btrig_ref(phases, refC, refS);

    auto report = [&](const char* label, auto eval){
        ErrAgg agg;
        for (int i=0;i<NUM_ANCHORS_FULL;++i) {
            long double c,s; eval(i, c, s);
            agg.add(c, s, refC[i], refS[i]);
        }
        agg.finalize(NUM_ANCHORS_FULL);
        std::cout << label
                  << "  rms|dcos|=" << agg.rmsC
                  << "  max|dcos|=" << agg.maxC
                  << "  rms|dsin|=" << agg.rmsS
                  << "  max|dsin|=" << agg.maxS << "\n";
    };

    report("btrig precise:", [&](int i, long double& c, long double& s){
        btrig::sincos_q64(phases[i], true, c, s);
    });
    report("btrig fast:   ", [&](int i, long double& c, long double& s){
        btrig::sincos_q64(phases[i], false, c, s);
    });
    report("sincos dbl:   ", [&](int i, long double& c, long double& s){
        double sd, cd; p_sincos((double)((long double)i / (long double)NUM_ANCHORS_FULL * btrig::TAU), &sd, &cd);
        c = (long double)cd; s = (long double)sd;
    });
    report("sincosl:      ", [&](int i, long double& c, long double& s){
        long double ang = ((long double)i / (long double)NUM_ANCHORS_FULL) * btrig::TAU;
        long double sd, cd; p_sincosl(ang, &sd, &cd);
        c = cd; s = sd;
    });
    std::cout << "\n";
}

// ------------------------------------------------------------------
// Boundary probes vs btrig precise (centered pivot crossings)
// ------------------------------------------------------------------
static void test_boundaries_vs_btrig() {
    std::cout << "=== Boundary probes vs btrig precise (centered pivot crossings) ===\n";
    const long double half_step_turns = 0.5L / (long double)NUM_ANCHORS_FULL;
    const long double one_step_turns  = 1.0L / (long double)NUM_ANCHORS_FULL;
    const long double eps = 1e-16L;

    std::vector<btrig::uq64> phases;
    phases.reserve(NUM_ANCHORS_FULL * 4);

    for (int i=0;i<NUM_ANCHORS_FULL;++i) {
        long double base = (long double)i / (long double)NUM_ANCHORS_FULL;
        long double pts[4] = {
            base,
            base + half_step_turns - eps,
            base + half_step_turns + eps,
            base + one_step_turns  - eps
        };
        for (long double ph : pts) {
            ph -= std::floor(ph);
            phases.push_back(turns_to_q64_rn(ph));
        }
    }

    std::vector<long double> refC, refS;
    compute_btrig_ref(phases, refC, refS);

    auto report = [&](const char* label, auto eval){
        ErrAgg agg;
        for (size_t i=0;i<phases.size();++i) {
            long double c,s; eval(i, c, s);
            agg.add(c, s, refC[i], refS[i]);
        }
        agg.finalize(phases.size());
        std::cout << label
                  << "  rms|dcos|=" << agg.rmsC
                  << "  max|dcos|=" << agg.maxC
                  << "  rms|dsin|=" << agg.rmsS
                  << "  max|dsin|=" << agg.maxS << "\n";
    };

    report("btrig precise:", [&](size_t i, long double& c, long double& s){
        btrig::sincos_q64(phases[i], true, c, s);
    });
    report("btrig fast:   ", [&](size_t i, long double& c, long double& s){
        btrig::sincos_q64(phases[i], false, c, s);
    });
    report("sincos dbl:   ", [&](size_t i, long double& c, long double& s){
        double sd, cd; p_sincos((double)q64_to_radians(phases[i]), &sd, &cd);
        c = (long double)cd; s = (long double)sd;
    });
    report("sincosl:      ", [&](size_t i, long double& c, long double& s){
        long double sd, cd; p_sincosl(q64_to_radians(phases[i]), &sd, &cd);
        c = cd; s = sd;
    });
    std::cout << "\n";
}

// ------------------------------------------------------------------
// Random phases: error vs btrig precise
// ------------------------------------------------------------------
static void test_random_vs_btrig(int N, uint64_t seed=0x5EEDBEEFULL) {
    auto phases_q = make_random_q64((size_t)N, seed);

    std::vector<long double> refC, refS;
    compute_btrig_ref(phases_q, refC, refS);

    auto eval_err = [&](auto eval){
        ErrAgg agg;
        for (int i=0;i<N;++i) {
            long double c,s; eval(i, c, s);
            agg.add(c, s, refC[i], refS[i]);
        }
        agg.finalize(N);
        return agg;
    };

    auto agg_precise = eval_err([&](int i, long double& c, long double& s){
        btrig::sincos_q64(phases_q[i], true, c, s);
    });
    auto agg_fast = eval_err([&](int i, long double& c, long double& s){
        btrig::sincos_q64(phases_q[i], false, c, s);
    });
    auto agg_sincos = eval_err([&](int i, long double& c, long double& s){
        double sd, cd; p_sincos((double)q64_to_radians(phases_q[i]), &sd, &cd);
        c = (long double)cd; s = (long double)sd;
    });
    auto agg_sincosl = eval_err([&](int i, long double& c, long double& s){
        long double sd, cd; p_sincosl(q64_to_radians(phases_q[i]), &sd, &cd);
        c = cd; s = sd;
    });

    std::cout << "=== Random phases: error vs btrig precise (N="<<N<<") ===\n"
              << std::setprecision(18) << std::fixed
              << "btrig precise: rms|dcos|=" << agg_precise.rmsC
              << "   max|dcos|=" << agg_precise.maxC
              << "   rms|dsin|=" << agg_precise.rmsS
              << "   max|dsin|=" << agg_precise.maxS << "\n"
              << "btrig fast:    rms|dcos|=" << agg_fast.rmsC
              << "   max|dcos|=" << agg_fast.maxC
              << "   rms|dsin|=" << agg_fast.rmsS
              << "   max|dsin|=" << agg_fast.maxS << "\n"
              << "sincos dbl:    rms|dcos|=" << agg_sincos.rmsC
              << "   max|dcos|=" << agg_sincos.maxC
              << "   rms|dsin|=" << agg_sincos.rmsS
              << "   max|dsin|=" << agg_sincos.maxS << "\n"
              << "sincosl:       rms|dcos|=" << agg_sincosl.rmsC
              << "   max|dcos|=" << agg_sincosl.maxC
              << "   rms|dsin|=" << agg_sincosl.rmsS
              << "   max|dsin|=" << agg_sincosl.maxS << "\n\n";
}

// ------------------------------------------------------------------
// Combined benches (timing) with errors vs btrig precise
// ------------------------------------------------------------------
static void bench_all(int N, uint64_t seed=0xBEEFFACEULL) {
    auto phases_q = make_random_q64((size_t)N, seed);

    std::vector<long double> refC, refS;
    compute_btrig_ref(phases_q, refC, refS);

    auto r1 = bench_btrig_q64<true >(phases_q, refC, refS);   // baseline
    auto r2 = bench_btrig_q64<false>(phases_q, refC, refS);   // vs baseline
    auto r3 = bench_sincos_vs_btrig (phases_q, refC, refS);   // vs baseline
    auto r4 = bench_sincosl_vs_btrig(phases_q, refC, refS);   // vs baseline

    auto print = [](const BenchResult& r){
        std::cout << r.name << "\n"
                  << "  time: " << std::setprecision(6) << std::fixed << r.secs << " s"
                  << "   ns/call: " << std::setprecision(2) << (r.ns_per_call)
                  << "   Mevals/s: " << std::setprecision(2) << r.mevals_per_s << "\n"
                  << "  RMS err cos/sin: " << std::setprecision(18) << r.rms_err_cos
                  << " / " << r.rms_err_sin << "\n"
                  << "  Max err cos/sin: " << r.max_err_cos
                  << " / " << r.max_err_sin << "\n"
                  << "  sink: " << r.sink << "\n\n";
    };

    print(r1); print(r2); print(r3); print(r4);
}

// ------------------------------------------------------------------

int main() {
    std::cout.setf(std::ios::fixed);
    std::cout << std::setprecision(18);

    // Accuracy suites: everything compared to btrig precise
    test_anchors_vs_btrig();
    test_random_vs_btrig(20000);
    test_boundaries_vs_btrig();

    // Speed benches (random phases), errors vs btrig precise
    bench_all(500000);

    // Touch the global sink so the linker cannot drop it
    std::cerr << "[sink guard] g_sink=" << (long double)g_sink << "\n";
    return 0;
}
