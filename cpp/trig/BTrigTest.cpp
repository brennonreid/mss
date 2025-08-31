// test_trig.cpp
#include <chrono>
#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include <vector>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <queue>
#include <limits>
#include <cstdio>
#include <cfloat>
#include <inttypes.h>  // for PRIu64

#include "btrig.h"
#include "btrigsimd.h"
#include "btrigdouble.h"   // <<< NEW: alternate accuracy path (hi/lo + optional DD)
#include "logexp.h"

using namespace std;
using namespace btrig;

// --- Best-of-3 control ---
static bool   g_quiet_prints = false;      // suppress prints for trial runs
static double g_best_ms = 1e300;           // best time across quiet trials
static volatile double g_best_sink = 0.0;  // matching sink for best time

static double print_bench(const char* name, double ms, volatile double sink) {
    if (g_quiet_prints) {
        if (ms < g_best_ms) {
            g_best_ms = ms;
            g_best_sink = sink;
        }
        return ms;
    }
    cout << left << setw(32) << name
        << " : " << fixed << setprecision(12) << setw(16) << ms
        << " ms    (sink=" << scientific << setprecision(17) << sink << ")\n"
        << defaultfloat;
    return ms;
}

// Per-call style (kept for continuity) but driven by arrays to unify inputs
template<typename F>
double bench_over_arrays_calls(const char* name, F func,
    const double* A, const double* B, int N) {
    volatile double sink = 0.0;
    auto t0 = chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i) sink += func(A[i], B[i]);
    auto t1 = chrono::high_resolution_clock::now();
    return print_bench(name,
        chrono::duration<double, milli>(t1 - t0).count(), sink);
}

// Array scalar baseline: one tight scalar loop calling btrig::sincos
double bench_array_scalar_sincos(const char* name, bool precise,
    const double* A, int N) {
    vector<double> C(N), S(N);
    volatile double sink = 0.0;
    auto t0 = chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i) btrig::sincos(A[i], precise, C[i], S[i]);
    auto t1 = chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i) sink += C[i] * 1.0000001 + S[i] * 0.9999999;
    return print_bench(name,
        chrono::duration<double, milli>(t1 - t0).count(), sink);
}

// SIMD batch (falls back to scalar if AVX2 disabled)
double bench_batch_sincos(const char* name, bool precise,
    const double* A, int N) {
    vector<double> C(N), S(N);
    volatile double sink = 0.0;
    auto t0 = chrono::high_resolution_clock::now();
    btrig::sincos_batch_avx2(A, N, precise, C.data(), S.data());
    auto t1 = chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i) sink += C[i] * 1.0000001 + S[i] * 0.9999999;
#if defined(__AVX2__)
    string label = string(name) + " (AVX2)";
#else
    string label = string(name) + " (scalar)";
#endif
    return print_bench(label.c_str(),
        chrono::duration<double, milli>(t1 - t0).count(), sink);
}

// --- Macro: best-of-3 without a fourth run ---
#define BEST_OF_3(_label_literal, _expr) do {         \
    g_best_ms = 1e300;                                 \
    g_best_sink = 0.0;                                 \
    g_quiet_prints = true;                             \
    for (int _k = 0; _k < 3; ++_k) { (void)(_expr); }  \
    g_quiet_prints = false;                            \
    print_bench((_label_literal), g_best_ms, g_best_sink); \
} while(0)

// ====== ULP helpers ======
static inline uint64_t ordered_bits(double x) {
    uint64_t u;
    std::memcpy(&u, &x, sizeof(double));
    return (u & 0x8000000000000000ull) ? ~u : (u ^ 0x8000000000000000ull);
}
static inline uint64_t ulp_diff(double a, double b) {
    if (!std::isfinite(a) || !std::isfinite(b)) return 0;
    uint64_t oa = ordered_bits(a), ob = ordered_bits(b);
    return (oa > ob) ? (oa - ob) : (ob - oa);
}
static inline double ulp_at(double x) {
    return std::nextafter(x, +INFINITY) - x;
}

// ====== ULP fairness gate ======
static inline bool ulp_is_fair(double refv, double /*gotv*/) {
    constexpr double ABS_ZERO_CUTOFF = 1e-13;
    const double ULP_AT_ONE = std::nextafter(1.0, +INFINITY) - 1.0;
    constexpr double MIN_ULP_FRACTION = 0.5;
    constexpr double K_ULP_NEAR_ZERO = 8.0;

    const double mag = std::fabs(refv);
    if (mag < ABS_ZERO_CUTOFF) return false;

    const double ulp_here = ulp_at(refv);
    if (ulp_here < ULP_AT_ONE * MIN_ULP_FRACTION) return false;
    if (mag < K_ULP_NEAR_ZERO * ulp_here) return false;
    return true;
}

// Optional: correctly-rounded refs if platform long double > 53 bits
static inline double ref_sin(double x) {
#if LDBL_MANT_DIG > 53
    return (double)sinl((long double)x);
#else
    return std::sin(x);
#endif
}
static inline double ref_cos(double x) {
#if LDBL_MANT_DIG > 53
    return (double)cosl((long double)x);
#else
    return std::cos(x);
#endif
}

// ====== Accuracy aggregation ======
struct ErrStats {
    double   max_abs = 0.0;
    double   rms_abs = 0.0;
    int      arg_index = -1;
    uint64_t n = 0;
    uint64_t sum_ulp = 0;
    uint64_t max_ulp = 0;
    double   avg_ulp = 0.0;
};

static inline void bump_ulp(ErrStats& es, uint64_t u) {
    ++es.n;
    es.sum_ulp += u;
    if (u > es.max_ulp) es.max_ulp = u;
}
static inline void finish_stats(ErrStats& es, long double sumsq, uint64_t denom) {
    es.rms_abs = std::sqrt((double)(sumsq / (long double)denom));
    es.avg_ulp = es.n ? (double)es.sum_ulp / (double)es.n : 0.0;
}
static inline void print_stats(const char* label, const ErrStats& es) {
    std::cout << std::left << std::setw(34) << label
        << " max_abs=" << std::scientific << std::setprecision(17) << es.max_abs
        << "  rms=" << std::scientific << std::setprecision(17) << es.rms_abs
        << "  max_ulp=" << std::defaultfloat << es.max_ulp
        << "  avg_ulp=" << std::fixed << std::setprecision(3) << es.avg_ulp
        << std::defaultfloat
        << "  N_ulp=" << es.n
        << "  @i=" << es.arg_index
        << "\n";
}

// sin/cos pair accuracy vs std:: (btrig baseline)
template<bool PRECISE>
ErrStats accuracy_sincos(const double* A, int N) {
    ErrStats es;
    long double sumsq = 0.0L;
    for (int i = 0; i < N; ++i) {
        double c_ref = std::cos(A[i]), s_ref = std::sin(A[i]);
        double c_got, s_got; btrig::sincos(A[i], PRECISE, c_got, s_got);

        double eC = std::fabs(c_got - c_ref);
        double eS = std::fabs(s_got - s_ref);
        sumsq += (long double)eC * (long double)eC;
        sumsq += (long double)eS * (long double)eS;

        double e = (eC > eS) ? eC : eS;
        if (e > es.max_abs) { es.max_abs = e; es.arg_index = i; }

        if (ulp_is_fair(c_ref, c_got)) bump_ulp(es, ulp_diff(c_got, c_ref));
        if (ulp_is_fair(s_ref, s_got)) bump_ulp(es, ulp_diff(s_got, s_ref));
    }
    finish_stats(es, sumsq, /*denom=*/2ull * (uint64_t)N);
    return es;
}

// ----- NEW: sin/cos pair accuracy vs std:: (btrigdouble) -----
template<bool PRECISE>
ErrStats accuracy_sincos_double(const double* A, int N) {
    ErrStats es;
    long double sumsq = 0.0L;
    for (int i = 0; i < N; ++i) {
        double c_ref = std::cos(A[i]), s_ref = std::sin(A[i]);
        double c_got, s_got; btrigdouble::sincos(A[i], PRECISE, c_got, s_got);

        double eC = std::fabs(c_got - c_ref);
        double eS = std::fabs(s_got - s_ref);
        sumsq += (long double)eC * (long double)eC;
        sumsq += (long double)eS * (long double)eS;

        double e = (eC > eS) ? eC : eS;
        if (e > es.max_abs) { es.max_abs = e; es.arg_index = i; }

        if (ulp_is_fair(c_ref, c_got)) bump_ulp(es, ulp_diff(c_got, c_ref));
        if (ulp_is_fair(s_ref, s_got)) bump_ulp(es, ulp_diff(s_got, s_ref));
    }
    finish_stats(es, sumsq, /*denom=*/2ull * (uint64_t)N);
    return es;
}

// single-output accuracy vs std::
template<typename FSTD, typename FBT>
ErrStats accuracy_single(const double* A, const double* B, int N, FSTD fstd, FBT fbt) {
    ErrStats es;
    long double sumsq = 0.0L;
    for (int i = 0; i < N; ++i) {
        double ref = fstd(A[i], B ? B[i] : 0.0);
        double got = fbt(A[i], B ? B[i] : 0.0);

        double e = std::fabs(got - ref);
        sumsq += (long double)e * (long double)e;
        if (e > es.max_abs) { es.max_abs = e; es.arg_index = i; }

        if (ulp_is_fair(ref, got)) bump_ulp(es, ulp_diff(got, ref));
    }
    finish_stats(es, sumsq, /*denom=*/(uint64_t)N);
    return es;
}

// ===== Random precise diffs =====
static void print_random_cs_diffs_precise(std::mt19937_64& rng, int count = 10) {
    std::uniform_real_distribution<double> dist(-TAU, TAU);

    std::cout << "\n--- RANDOM ANGLE CHECK (std vs btrig precise, " << count << " samples) ---\n";
    std::cout << std::left
        << std::setw(3) << "#"
        << std::setw(24) << "angle (rad)"
        << std::setw(24) << "std cos"
        << std::setw(24) << "std sin"
        << std::setw(24) << "btrig cos (prec)"
        << std::setw(24) << "btrig sin (prec)"
        << std::setw(24) << "|cos diff|"
        << std::setw(24) << "|sin diff|"
        << "\n";

    for (int i = 0; i < count; ++i) {
        double a = dist(rng);
        double c_ref = std::cos(a);
        double s_ref = std::sin(a);
        double c_prec, s_prec;
        btrig::sincos(a, /*precise=*/true, c_prec, s_prec);

        std::cout << std::right << std::setw(3) << (i + 1)
            << std::fixed << std::setprecision(17)
            << std::setw(24) << a
            << std::setw(24) << c_ref
            << std::setw(24) << s_ref
            << std::setw(24) << c_prec
            << std::setw(24) << s_prec
            << std::setw(24) << std::fabs(c_prec - c_ref)
            << std::setw(24) << std::fabs(s_prec - s_ref)
            << "\n";
    }
}

enum class StepSource { BTRIG, STD };
static inline void get_step_pair(double delta, bool precise, StepSource src, double& cd, double& sd) {
    if (src == StepSource::BTRIG) {
        btrig::sincos(delta, precise, cd, sd);
    }
    else {
        cd = std::cos(delta);
        sd = std::sin(delta);
    }
}

// ===== Audits: floor/residual logic and LUT continuity/unit length =====
namespace btrig {

    static inline int floor_idx_ref(double t) {
        return (int)std::floor(t);
    }
    static inline int floor_idx_fast(double t) {
        int i = (int)t;
        i -= ((t < 0.0) & ((double)i != t));
        return i;
    }

    static void audit_floor_and_residual() {
        std::puts("=== AUDIT: floor/residual ===");
        auto chk = [&](double ang) {
            const double t = ang * ONE_OVER_STEP;
            const int i_ref = floor_idx_ref(t);
            const int i_fast = floor_idx_fast(t);

            if (i_ref != i_fast) {
                std::printf("FLOOR MISMATCH  t=%.17g  i_ref=%d  i_fast=%d\n", t, i_ref, i_fast);
            }

            const double d_ref = std::fma(-(double)i_ref, STEP, ang);
            const double d_fast = std::fma(-(double)i_fast, STEP, ang);
            if (d_ref != d_fast) {
                std::printf("RESIDUAL MISMATCH ang=%.17g  d_ref=%.17g  d_fast=%.17g\n",
                    ang, d_ref, d_fast);
            }
            };
        for (int k = -10; k <= 10; ++k) {
            double ang = (115 + 0.5 + 0.02 * k) * STEP;
            chk(ang);
        }
        for (int m = 0; m <= 6; ++m) {
            double ang = -(double)m * STEP;
            chk(ang);
        }
        std::puts("=== END AUDIT: floor/residual ===");
    }

    static void audit_lut() {
        std::puts("=== AUDIT: LUT radius and continuity ===");
        const int N = NUM_ANCHORS_FULL;
        const double ideal_dp = std::cos(STEP);
        auto dot = [](double c0, double s0, double c1, double s1) {
            return c0 * c1 + s0 * s1;
            };

        int bad_len = 0, bad_dp = 0;
        for (int k = 0; k < N; ++k) {
            const double c = fullCircleAnchors[k][1];
            const double s = fullCircleAnchors[k][2];
            const double r2 = std::fma(c, c, s * s);
            const double len_err = std::abs(r2 - 1.0);
            if (len_err > 1e-15) {
                ++bad_len;
                std::printf("LEN != 1  idx=%d  err=%.3g  c=%.17g  s=%.17g\n", k, len_err, c, s);
            }
            const int kn = (k + 1) & (N - 1);
            const double cn = fullCircleAnchors[kn][1];
            const double sn = fullCircleAnchors[kn][2];
            const double dp = dot(c, s, cn, sn);
            const double seam_err = std::abs(dp - ideal_dp);
            if (seam_err > 1e-15 && (k == 115 || k == 116)) {
                ++bad_dp;
                std::printf("SEAM dp mismatch at %d->%d  dp-ideal=%.3g  dp=%.17g  ideal=%.17g\n",
                    k, kn, seam_err, dp, ideal_dp);
            }
        }
        if (!bad_len) std::puts("LUT radius OK (all anchors unit within 1e-15)");
        if (!bad_dp)  std::puts("Seam neighbor dot OK within 1e-15 at 115/116");
        std::puts("=== END AUDIT: LUT ===");
    }

    static void dump_seam_window(int center_idx = 115, int halfspan = 2) {
        std::puts("=== DUMP: seam window ===");
        for (int dk = -halfspan; dk <= halfspan; ++dk) {
            int k = (center_idx + dk) & (NUM_ANCHORS_FULL - 1);
            double c = fullCircleAnchors[k][1];
            double s = fullCircleAnchors[k][2];
            std::printf("idx=%3d  c=%.17g  s=%.17g  r2=%.17g\n",
                k, c, s, std::fma(c, c, s * s));
        }
        std::puts("=== END DUMP ===");
    }
    static void run_audits() {
        audit_floor_and_residual();
        audit_lut();
        dump_seam_window(115, 3);
    }
} // namespace btrig

// ====== ULP trouble-maker capture ======
struct PathInfo { int idx = 0; double d = 0.0; double a1 = 0.0; double a2 = 0.0; };
static inline PathInfo decode_edge(double angle) {
    PathInfo P{};
    double t = angle * btrig::ONE_OVER_STEP;
    int i = (int)std::floor(t);
    P.idx = i & btrig::ANCHOR_MASK;
    P.d = (t - i) * btrig::STEP;
    P.a1 = btrig::fullCircleAnchors[P.idx][1];
    P.a2 = btrig::fullCircleAnchors[P.idx][2];
    return P;
}
static inline PathInfo decode_center(double angle) {
    PathInfo P{};
    double t = angle * btrig::ONE_OVER_STEP;
    int i = (int)(t + (t >= 0.0 ? 0.5 : -0.5));
    if (t < 0.0 && (double)i > t) --i;
    P.idx = i & btrig::ANCHOR_MASK;
    P.d = (t - i) * btrig::STEP;
    P.a1 = btrig::fullCircleAnchors[P.idx][1];
    P.a2 = btrig::fullCircleAnchors[P.idx][2];
    return P;
}
struct Hit {
    bool is_sin = false;
    int i = -1;
    double angle = 0.0;
    uint64_t ulp = 0;
    double abs_err = 0.0;
    double got = 0.0;
    double refv = 0.0;
    double ulp_size = 0.0;
    PathInfo edge;
    PathInfo center;
};
static std::vector<Hit> topk_trouble_sincos(const double* A, int N, int K) {
    struct Cmp { bool operator()(const Hit& a, const Hit& b) const { return a.ulp > b.ulp; } };
    std::priority_queue<Hit, std::vector<Hit>, Cmp> pq;
    auto consider = [&](bool is_sin, int i, double angle, double got, double refv) {
        if (!ulp_is_fair(refv, got)) return;
        uint64_t u = ulp_diff(got, refv);
        Hit h{}; h.is_sin = is_sin; h.i = i; h.angle = angle; h.ulp = u;
        h.abs_err = std::fabs(got - refv);
        h.got = got; h.refv = refv;
        h.ulp_size = ulp_at(refv);
        h.edge = decode_edge(angle);
        h.center = decode_center(angle);
        if ((int)pq.size() < K) pq.push(h);
        else if (u > pq.top().ulp) { pq.pop(); pq.push(h); }
        };
    for (int i = 0; i < N; ++i) {
        double a = A[i];
        double c_ref = ref_cos(a), s_ref = ref_sin(a);
        double c_got, s_got; btrig::sincos(a, /*precise=*/true, c_got, s_got);
        consider(false, i, a, c_got, c_ref);
        consider(true, i, a, s_got, s_ref);
    }
    std::vector<Hit> out; out.reserve(pq.size());
    while (!pq.empty()) { out.push_back(pq.top()); pq.pop(); }
    std::sort(out.begin(), out.end(), [](const Hit& A, const Hit& B) { return A.ulp > B.ulp; });
    return out;
}

// ----- NEW: TOPK trouble-maker for btrigdouble -----
static std::vector<Hit> topk_trouble_sincos_double(const double* A, int N, int K) {
    struct Cmp { bool operator()(const Hit& a, const Hit& b) const { return a.ulp > b.ulp; } };
    std::priority_queue<Hit, std::vector<Hit>, Cmp> pq;
    auto consider = [&](bool is_sin, int i, double angle, double got, double refv) {
        if (!ulp_is_fair(refv, got)) return;
        uint64_t u = ulp_diff(got, refv);
        Hit h{}; h.is_sin = is_sin; h.i = i; h.angle = angle; h.ulp = u;
        h.abs_err = std::fabs(got - refv);
        h.got = got; h.refv = refv;
        h.ulp_size = ulp_at(refv);
        h.edge = decode_edge(angle);
        h.center = decode_center(angle);
        if ((int)pq.size() < K) pq.push(h);
        else if (u > pq.top().ulp) { pq.pop(); pq.push(h); }
        };
    for (int i = 0; i < N; ++i) {
        double a = A[i];
        double c_ref = ref_cos(a), s_ref = ref_sin(a);
        double c_got, s_got; btrigdouble::sincos(a, /*precise=*/true, c_got, s_got);
        consider(false, i, a, c_got, c_ref);
        consider(true, i, a, s_got, s_ref);
    }
    std::vector<Hit> out; out.reserve(pq.size());
    while (!pq.empty()) { out.push_back(pq.top()); pq.pop(); }
    std::sort(out.begin(), out.end(), [](const Hit& A, const Hit& B) { return A.ulp > B.ulp; });
    return out;
}

static void print_trouble_list(const std::vector<Hit>& hits) {
    std::cout << "\n=== TOP TROUBLE CASES (sorted by ULP error) ===\n";
    std::cout << std::left
        << std::setw(5) << "#"
        << std::setw(6) << "func"
        << std::setw(23) << "angle_rad"
        << std::setw(14) << "ulp_err"
        << std::setw(22) << "abs_err"
        << std::setw(23) << "ref_val"
        << std::setw(23) << "got_val"
        << std::setw(7) << "L_idx"
        << std::setw(14) << "L_d(rad)"
        << std::setw(12) << "L_a1"
        << std::setw(12) << "L_a2"
        << std::setw(7) << "C_idx"
        << std::setw(14) << "C_d(rad)"
        << std::setw(12) << "C_a1"
        << std::setw(12) << "C_a2"
        << std::setw(14) << "ulp_unit@ref"
        << "\n";
    std::cout << std::string(5 + 6 + 23 + 14 + 22 + 23 + 23 + 7 + 14 + 12 + 12 + 7 + 14 + 12 + 12 + 14, '-') << "\n";
    for (size_t k = 0; k < hits.size(); ++k) {
        const Hit& h = hits[k];
        std::cout << std::right
            << std::setw(5) << (k + 1)
            << std::setw(6) << (h.is_sin ? "sin" : "cos")
            << std::fixed << std::setprecision(17)
            << std::setw(23) << h.angle
            << std::defaultfloat
            << std::setw(14) << h.ulp
            << std::scientific << std::setprecision(6)
            << std::setw(22) << h.abs_err
            << std::fixed << std::setprecision(17)
            << std::setw(23) << h.refv
            << std::setw(23) << h.got
            << std::defaultfloat
            << std::setw(7) << h.edge.idx
            << std::scientific << std::setprecision(6)
            << std::setw(14) << h.edge.d
            << std::fixed << std::setprecision(9)
            << std::setw(12) << h.edge.a1
            << std::setw(12) << h.edge.a2
            << std::defaultfloat
            << std::setw(7) << h.center.idx
            << std::scientific << std::setprecision(6)
            << std::setw(14) << h.center.d
            << std::fixed << std::setprecision(9)
            << std::setw(12) << h.center.a1
            << std::setw(12) << h.center.a2
            << std::scientific << std::setprecision(6)
            << std::setw(14) << h.ulp_size
            << std::defaultfloat
            << "\n";
    }
}

// ====== Tracer (current path vs hardened control) ======
static inline int fast_floor_idx(double t) {
    int i = (int)t;
    i -= ((t < 0.0) & ((double)i != t));
    return i;
}
static inline void hardened_index_residual(double angle, int& idx, double& d) {
    const double k = angle / TAU;
    long long q = (long long)k;
    if (angle < 0.0 && (double)q != k) --q;

    const double ang = std::fma(-(double)q, TAU, angle);

    const double t = ang * ONE_OVER_STEP;
    int i = fast_floor_idx(t);
    idx = i & ANCHOR_MASK;

    d = std::fma(-(double)i, STEP, ang);

    // wrap using the mask (since NUM_ANCHORS_FULL is a power of two)
    if (d >= STEP) { d -= STEP; idx = (idx + 1) & ANCHOR_MASK; }
    else if (d < 0.0) { d += STEP; idx = (idx - 1) & ANCHOR_MASK; }
}

static inline double len2(double x, double y) { return std::fma(x, x, y * y); }

static void trace_angle(double angle, bool precise) {
    using std::printf;
    puts("\n=== TRACE ANGLE ===");
    printf("angle = %.17g\n", angle);

    // current path
    const double t = angle * ONE_OVER_STEP;
    const int i_fast = fast_floor_idx(t);
    const int idx_fast = i_fast & ANCHOR_MASK;
    const double d_fast = std::fma(-(double)i_fast, STEP, angle);

    const double* aC = fullCircleAnchors[idx_fast];
    const double cA = aC[1], sA = aC[2];

    const double d2 = d_fast * d_fast;
    double dx = 0.0, dy = 0.0;
    if (precise) {
        double p = std::fma(-d2, INV_40320, INV_720);
        p = std::fma(-d2, p, INV_24);
        p = std::fma(-d2, p, INV_2);
        dx = std::fma(-d2, p, 1.0);
        double q = std::fma(-d2, INV_5040, INV_120);
        q = std::fma(-d2, q, INV_6);
        dy = d_fast * std::fma(-d2, q, 1.0);
    }
    else {
        dx = std::fma(-d2, INV_2, 1.0);
        dy = d_fast;
    }
    double cosB = std::fma(dx, cA, -dy * sA);
    double sinB = std::fma(dx, sA, dy * cA);

    const double cosS = std::cos(angle);
    const double sinS = std::sin(angle);

    const uint64_t ulpC = ulp_diff(cosB, cosS);
    const uint64_t ulpS = ulp_diff(sinB, sinS);
    printf("[current] t=%.17g  i=%d  idx=%d  d=%.17g\n", t, i_fast, idx_fast, d_fast);
    printf("[anchor ] c=%.17g  s=%.17g  r2=%.17g\n", cA, sA, len2(cA, sA));
    printf("[micro  ] dx=%.17g  dy=%.17g  len2=%.17g\n", dx, dy, len2(dx, dy));
    printf("[result ] cos=%.17g  sin=%.17g\n", cosB, sinB);
    printf("[std    ] cos=%.17g  sin=%.17g\n", cosS, sinS);
    printf("[ULP    ] cos=%" PRIu64 "  sin=%" PRIu64 "    |dc|=%.3e |ds|=%.3e\n",
        ulpC, ulpS, std::fabs(cosB - cosS), std::fabs(sinB - sinS));

    // control path
    int idxH = 0; double dH = 0.0; hardened_index_residual(angle, idxH, dH);
    const double* aH = fullCircleAnchors[idxH];
    const double d2H = dH * dH;
    double dxH = 0.0, dyH = 0.0;
    if (precise) {
        double p = std::fma(-d2H, INV_40320, INV_720);
        p = std::fma(-d2H, p, INV_24);
        p = std::fma(-d2H, p, INV_2);
        dxH = std::fma(-d2H, p, 1.0);
        double q = std::fma(-d2H, INV_5040, INV_120);
        q = std::fma(-d2H, q, INV_6);
        dyH = dH * std::fma(-d2H, q, 1.0);
    }
    else {
        dxH = std::fma(-d2H, INV_2, 1.0);
        dyH = dH;
    }
    const double cosH = std::fma(dxH, aH[1], -dyH * aH[2]);
    const double sinH = std::fma(dxH, aH[2], dyH * aH[1]);
    const uint64_t ulpCH = ulp_diff(cosH, cosS);
    const uint64_t ulpSH = ulp_diff(sinH, sinS);
    printf("[control] idx=%d  d=%.17g  (c,s)=(%.17g,%.17g)\n", idxH, dH, aH[1], aH[2]);
    printf("[control] cos=%.17g  sin=%.17g  ULPc=%" PRIu64 "  ULPs=%" PRIu64 "\n",
        cosH, sinH, ulpCH, ulpSH);
    printf("[delta   ] dcos=%.3e  dsin=%.3e  (current - control)\n",
        cosB - cosH, sinB - sinH);

    // small sweep around the angle to visualize seam behavior
    const double step = STEP * 0.02;
    puts("\n[mini-sweep around angle]");
    puts("   k    d(rad)            ULPc   ULPs");
    for (int k = -10; k <= 10; ++k) {
        double a2 = angle + k * step;
        double cG, sG; btrig::sincos(a2, precise, cG, sG);
        uint64_t uc = ulp_diff(cG, std::cos(a2));
        uint64_t us = ulp_diff(sG, std::sin(a2));
        printf("%4d  %+ .17g   %6" PRIu64 "  %6" PRIu64 "\n", k, k * step, uc, us);
    }
    puts("=== END TRACE ===");
}

// ====== MAIN ======
int main() {
#if defined(__AVX2__)
    std::cout << "AVX2 is enabled.\n";
#elif defined(__AVX__)
    std::cout << "AVX is enabled.\n";
#elif defined(__SSE4_1__)
    std::cout << "SSE4.1 is enabled.\n";
#elif defined(__SSE2__)
    std::cout << "SSE2 is enabled.\n";
#else
    std::cout << "No SIMD detected.\n";
#endif

    const int N = 15'500'000;

    // Build shared inputs ONCE so every test uses the same ranges/values
    mt19937_64 rng(12345);
    uniform_real_distribution<double> dist(-TAU, TAU);

    vector<double> A(N), B(N), Z(N);
    for (int i = 0; i < N; ++i) {
        A[i] = dist(rng);
        B[i] = dist(rng);
        Z[i] = B[i];
    }

    // NEW: initialize hi/lo pivots for btrigdouble once

    btrig::initFullCircleAnchors();
    btrigdouble::init();


    cout << "\n--- TRIG SPEED TEST (" << N << " samples, unified inputs) ---\n\n";

    BEST_OF_3("std::sin/cos pair (calls)",
        bench_over_arrays_calls("std::sin/cos pair (calls)", [](double a, double) {
            return std::sin(a) + std::cos(a);
            }, A.data(), B.data(), N));

    BEST_OF_3("btrigdouble::sincos (calls)",
        bench_over_arrays_calls("btrigdouble::sincos (calls)", [](double a, double) {
            double c, s; btrigdouble::sincos(a, true, c, s); return c + s;
            }, A.data(), B.data(), N));

    BEST_OF_3("btrig::sincos precise (calls)",
        bench_over_arrays_calls("btrig::sincos precise (calls)", [](double a, double) {
            double c, s; btrig::sincos(a, true, c, s); return c + s;
            }, A.data(), B.data(), N));

    // ----- NEW: btrigdouble accuracy path benches -----
    /*
    BEST_OF_3("btrigdouble::sincos fast (calls)",
        bench_over_arrays_calls("btrigdouble::sincos fast (calls)", [](double a, double) {
            double c, s; btrigdouble::sincos(a, false, c, s); return c + s;
            }, A.data(), B.data(), N));
    */
    

    BEST_OF_3("btrig::sincos fast (calls)",
        bench_over_arrays_calls("btrig::sincos fast (calls)", [](double a, double) {
            double c, s; btrig::sincos(a, false, c, s); return c + s;
            }, A.data(), B.data(), N));


    cout << "\n";

    BEST_OF_3("std::sin only (calls)",
        bench_over_arrays_calls("std::sin only (calls)", [](double a, double) {
            return std::sin(a);
            }, A.data(), B.data(), N));

    BEST_OF_3("btrig::sin fast (calls)",
        bench_over_arrays_calls("btrig::sin fast (calls)", [](double a, double) {
            double s; btrig::sin(a, false, s); return s;
            }, A.data(), B.data(), N));

    BEST_OF_3("btrig::sin precise (calls)",
        bench_over_arrays_calls("btrig::sin precise (calls)", [](double a, double) {
            double s; btrig::sin(a, true, s); return s;
            }, A.data(), B.data(), N));

    cout << "\n";

    BEST_OF_3("std::cos only (calls)",
        bench_over_arrays_calls("std::cos only (calls)", [](double a, double) {
            return std::cos(a);
            }, A.data(), B.data(), N));

    BEST_OF_3("btrig::cos fast (calls)",
        bench_over_arrays_calls("btrig::cos fast (calls)", [](double a, double) {
            double c; btrig::cos(a, false, c); return c;
            }, A.data(), B.data(), N));

    BEST_OF_3("btrig::cos precise (calls)",
        bench_over_arrays_calls("btrig::cos precise (calls)", [](double a, double) {
            double c; btrig::cos(a, true, c);  return c;
            }, A.data(), B.data(), N));

    cout << "\n";

    BEST_OF_3("btrig::sincos array fast (scalar)",
        bench_array_scalar_sincos("btrig::sincos array fast (scalar)", false, A.data(), N));

    BEST_OF_3("btrig::sincos array precise (scalar)",
        bench_array_scalar_sincos("btrig::sincos array precise (scalar)", true, A.data(), N));

#if defined(__AVX2__)
    const char* LBL_BATCH_FAST = "btrig::sincos_batch fast (AVX2)";
    const char* LBL_BATCH_PRECISE = "btrig::sincos_batch precise (AVX2)";
#else
    const char* LBL_BATCH_FAST = "btrig::sincos_batch fast (scalar)";
    const char* LBL_BATCH_PRECISE = "btrig::sincos_batch precise (scalar)";
#endif

    BEST_OF_3(LBL_BATCH_FAST,
        bench_batch_sincos("btrig::sincos_batch fast", false, A.data(), N));
    BEST_OF_3(LBL_BATCH_PRECISE,
        bench_batch_sincos("btrig::sincos_batch precise", true, A.data(), N));

    cout << "\n--- END TEST ---\n";

    // =========================
    // ACCURACY REPORT vs std::
    // =========================
    cout << "\n=== ACCURACY REPORT vs std:: ===\n";
    {
        //auto e_fast_d = accuracy_sincos_double<false>(A.data(), N);
        auto e_prec_d = accuracy_sincos_double<true >(A.data(), N);
        //print_stats("btrigdouble::sincos fast", e_fast_d);
        print_stats("btrigdouble::sincos precise", e_prec_d);
    }
    
    {
        auto e_fast = accuracy_sincos<false>(A.data(), N);
        auto e_precise = accuracy_sincos<true >(A.data(), N);
        print_stats("btrig::sincos precise", e_precise);
        print_stats("btrig::sincos fast", e_fast);
    }
    
    {
        auto e_sin_fast = accuracy_single(A.data(), nullptr, N,
            [](double a, double) { return std::sin(a); },
            [](double a, double) { double s; btrig::sin(a, false, s); return s; });
        auto e_sin_prec = accuracy_single(A.data(), nullptr, N,
            [](double a, double) { return std::sin(a); },
            [](double a, double) { double s; btrig::sin(a, true, s);  return s; });
        print_stats("btrig::sin precise", e_sin_prec);
        print_stats("btrig::sin fast", e_sin_fast);
        
    }
    {
        auto e_cos_fast = accuracy_single(A.data(), nullptr, N,
            [](double a, double) { return std::cos(a); },
            [](double a, double) { double c; btrig::cos(a, false, c); return c; });
        auto e_cos_prec = accuracy_single(A.data(), nullptr, N,
            [](double a, double) { return std::cos(a); },
            [](double a, double) { double c; btrig::cos(a, true, c);  return c; });
        print_stats("btrig::cos precise", e_cos_prec);
        print_stats("btrig::cos fast", e_cos_fast);
        
    }

    

    cout << "=== END ACCURACY REPORT ===\n";

    // Random sample dump (original btrig precise)
    print_random_cs_diffs_precise(rng, 10);

    // Audits (index/residual + LUT seam) for original btrig
    btrig::run_audits();

    // =========================
    // ULP TROUBLE-MAKER REPORTS + auto-trace
    // =========================
    {
        const int TOPK = 20;
        // Original
        auto worst = topk_trouble_sincos(A.data(), N, TOPK);
        print_trouble_list(worst);

        if (!worst.empty()) {
            const auto& h = worst.front();
            std::cout << "\n>>> Auto-tracing worst case (btrig "
                << (h.is_sin ? "sin" : "cos")
                << "), angle=" << std::setprecision(17) << h.angle
                << ", ULP=" << h.ulp << "\n";
            trace_angle(h.angle, /*precise=*/true);
        }

        cout << "=== btrigdouble ===\n";
        // NEW: btrigdouble
        auto worst_double = topk_trouble_sincos_double(A.data(), N, TOPK);
        print_trouble_list(worst_double);
    }

    return 0;
}
