#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include "TaylorTrig.h"

// Prevent the compiler from optimizing the loops away
static volatile double sink_double = 0.0;

static inline double bench_std_once(int N) {
    double checksum = 0.0;
    const double da = TAU / static_cast<double>(N);
    double a = 0.0;

    const auto t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i) {
        checksum += std::sin(a) + std::cos(a);
        a += da; // no manual wrap
    }
    const auto t1 = std::chrono::high_resolution_clock::now();

    sink_double = checksum;
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

static inline double bench_custom_pair_once(int N, bool precise) {
    double checksum = 0.0;
    const double da = TAU / static_cast<double>(N);
    double a = 0.0;

    const auto t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i) {
        double c, s;
        getUnitVectorFromAngle2(a, precise, c, s);
        checksum += c + s;
        a += da; // no manual wrap
    }
    const auto t1 = std::chrono::high_resolution_clock::now();

    sink_double = checksum;
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

static inline double bench_custom_singletons_once(int N, bool precise) {
    double checksum = 0.0;
    const double da = TAU / static_cast<double>(N);
    double a = 0.0;

    const auto t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i) {
        double s, c;
        customSin(a, precise, s);
        customCos(a, precise, c);
        checksum += c + s;
        a += da; // no manual wrap
    }
    const auto t1 = std::chrono::high_resolution_clock::now();

    sink_double = checksum;
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

// Repeat a few times and take the best (min) to reduce noise
static inline double repeat_min(double (*fn)(int, bool), int N, bool arg, int reps, double* out_checksum) {
    double best = 1e300;
    double best_checksum = 0.0;
    for (int i = 0; i < reps; ++i) {
        const double t = fn(N, arg);
        if (t < best) {
            best = t;
            best_checksum = sink_double; // grab checksum from THIS run
        }
    }
    if (out_checksum) *out_checksum = best_checksum;
    return best;
}

static inline double repeat_min_std(double (*fn)(int), int N, int reps, double* out_checksum) {
    double best = 1e300;
    double best_checksum = 0.0;
    for (int i = 0; i < reps; ++i) {
        const double t = fn(N);
        if (t < best) {
            best = t;
            best_checksum = sink_double;
        }
    }
    if (out_checksum) *out_checksum = best_checksum;
    return best;
}


int main() {
    const int N = 10000000;
    const int reps = 3;

    std::cout << "--- SPEED TEST (" << N << " samples, min of " << reps << ") ---\n";

    double chk1, chk2, chk3, chk4, chk5;
    const double t_std = repeat_min_std(bench_std_once, N, reps, &chk1);
    const double t_fast_pair = repeat_min(bench_custom_pair_once, N, false, reps, &chk2);
    const double t_prec_pair = repeat_min(bench_custom_pair_once, N, true, reps, &chk3);
    const double t_fast_sing = repeat_min(bench_custom_singletons_once, N, false, reps, &chk4);
    const double t_prec_sing = repeat_min(bench_custom_singletons_once, N, true, reps, &chk5);

    std::cout << "checksum1: " << chk1 << "\n";
    std::cout << "checksum2: " << chk2 << "\n";
    std::cout << "checksum3: " << chk3 << "\n";
    std::cout << "checksum4: " << chk4 << "\n";
    std::cout << "checksum5: " << chk5 << "\n\n";

    std::cout << std::setprecision(17) << "checksum1: " << chk1 << "\n";
    std::cout << std::setprecision(17) << "checksum3: " << chk3 << "\n";
    std::cout << std::setprecision(17) << "checksum5: " << chk5 << "\n";
    


    std::cout << std::setprecision(6) << "std::sin/cos (ms):            " << t_std << "\n";
    std::cout << std::setprecision(6) << "custom pair fast  (ms):       " << t_fast_pair << "  (getUnitVectorFromAngle2)\n";
    std::cout << std::setprecision(6) << "custom pair precise(ms):       " << t_prec_pair << "  (getUnitVectorFromAngle2)\n";
    std::cout << std::setprecision(6) << "custom singletons fast (ms):  " << t_fast_sing << "  (customSin + customCos)\n";
    std::cout << std::setprecision(6) << "custom singletons precise(ms): " << t_prec_sing << "  (customSin + customCos)\n\n";

    std::cout << "--- ACCURACY TEST (10 samples) ---\n";
    for (int i = 0; i <= 10; ++i) {
        const double a = (TAU * i) / 10.0;
        const double s_ref = std::sin(a), c_ref = std::cos(a);

        double cf, sf, cp, sp;
        getUnitVectorFromAngle2(a, false, cf, sf);
        getUnitVectorFromAngle2(a, true, cp, sp);

        std::cout << "angle=" << a
            << "  fast:    sin_err=" << (sf - s_ref)
            << "  cos_err=" << (cf - c_ref)
            << "  precise: sin_err=" << (sp - s_ref)
            << "  cos_err=" << (cp - c_ref)
            << "\n";
    }
    return 0;
}
