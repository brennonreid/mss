// BTrig.cpp : Seam scans on known difficult sin/cos locations (cardinals + diagonals)
// and a throughput benchmark against std::sin/std::cos

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <cfloat>

#include "btrig.h"

// -----------------------------
// Global timing clock alias
// -----------------------------
using Clock = std::chrono::high_resolution_clock;

// -----------------------------
// Seam-scan configuration
// -----------------------------

// Known difficult seams for sin/cos:
// - Cardinals:  0.00, 0.25, 0.50, 0.75  (0°,  90°, 180°, 270°)
// - Diagonals:  0.125, 0.375, 0.625, 0.875 (45°, 135°, 225°, 315°)
struct Seam {
    const char* label;  // short token for filenames/prints
    double      phase;  // turns
};

static const Seam SEAMS[] = {
    {"cardinal_0deg",    0.00000000000000000}, // 0°
    {"diagonal_45deg",   0.12500000000000000}, // 45°
    {"cardinal_90deg",   0.25000000000000000}, // 90°
    {"diagonal_135deg",  0.37500000000000000}, // 135°
    {"cardinal_180deg",  0.50000000000000000}, // 180°
    {"diagonal_225deg",  0.62500000000000000}, // 225°
    {"cardinal_270deg",  0.75000000000000000}, // 270°
    {"diagonal_315deg",  0.87500000000000000}, // 315°
};
static const int NUM_SEAMS = (int)(sizeof(SEAMS) / sizeof(SEAMS[0]));

// Neighborhood scan parameters
static const int    STEPS_PER_SIDE = 1000000;   // 1,000,000 either side
static const bool   USE_PRECISE = true;     // precise polynomial path
static const double PHASE_STEP = 1e-12;     // turns (tight neighborhood)

// Histogram + worst sample capture
static const int    HIST_BINS = 200;
static const int    TOP_KEEP = 1000;     // top N to write
static const int    RESERVOIR_SIZE = 2000;     // internal buffer for selecting top

// -----------------------------
// Throughput benchmark configuration
// -----------------------------
static const size_t NUM_SAMPLES = 15500000;    // ~15.5M like prior tests
static const double PHASE_SPAN = 1.0;         // sweep one full turn
static const double PHASE_OFFSET = 0.123456789; // de-align from cardinals
static const int    REPEATS = 1;           // increase to average multiple runs

// -----------------------------
// Helpers
// -----------------------------

// Build a deterministic, wrap-safe phase sequence in [0,1)
static std::vector<double> make_phases(size_t n) {
    std::vector<double> v; v.reserve(n);
    const double step = PHASE_SPAN / static_cast<double>(n);
    double p = PHASE_OFFSET;
    for (size_t i = 0; i < n; ++i) {
        p += step;
        // wrap to [0,1) without branching
        p -= std::floor(p);
        v.push_back(p);
    }
    return v;
}

static inline double now_ms() {
    return std::chrono::duration<double, std::milli>(Clock::now().time_since_epoch()).count();
}

struct SampleErr {
    double phase;
    double angle;      // radians
    double cos_b;      // btrig
    double cos_r;      // std
    double sin_b;      // btrig
    double sin_r;      // std
    double err_c_abs;
    double err_s_abs;
    double err_mag;    // sqrt(ec^2 + es^2)
};

static void scan_seam(const Seam& seam) {
    auto t0 = Clock::now();

    const long long total = 2LL * STEPS_PER_SIDE + 1LL;

    // Pretty banner
    const double center_phase = seam.phase;
    const double center_angle = center_phase * btrig::TAU;
    const double center_deg = center_phase * 360.0;

    std::cout.setf(std::ios::fixed); std::cout << std::setprecision(17);
    std::cout << "\n=== Seam: " << seam.label
        << " | phase=" << center_phase
        << " | angle_rad=" << center_angle
        << " | angle_deg=" << center_deg << " ===\n";

    // Accumulators
    double sum_sq_cos = 0.0, sum_sq_sin = 0.0;
    double max_cos_err = 0.0, max_cos_err_phase = 0.0;
    double max_sin_err = 0.0, max_sin_err_phase = 0.0;
    double global_max_mag = 0.0;

    std::vector<SampleErr> reservoir;
    reservoir.reserve(RESERVOIR_SIZE);

    // pass 1: stats, global max magnitude, worst reservoir selection
    for (long long k = -STEPS_PER_SIDE; k <= STEPS_PER_SIDE; ++k) {
        const double dk = static_cast<double>(k);
        const double phase = std::fma(dk, PHASE_STEP, center_phase);    // CHANGE: no drift
        const double angle = std::fma(phase, btrig::TAU, 0.0);          // CHANGE: fused scale

        double c_b, s_b;
        btrig::sincos(phase, USE_PRECISE, c_b, s_b);

        const double c_r = std::cos(angle);
        const double s_r = std::sin(angle);

        const double e_c = std::fabs(c_b - c_r);
        const double e_s = std::fabs(s_b - s_r);
        sum_sq_cos += e_c * e_c;
        sum_sq_sin += e_s * e_s;

        if (e_c > max_cos_err) { max_cos_err = e_c; max_cos_err_phase = phase; }
        if (e_s > max_sin_err) { max_sin_err = e_s; max_sin_err_phase = phase; }

        const double mag = std::sqrt(e_c * e_c + e_s * e_s);
        if (mag > global_max_mag) global_max_mag = mag;

        if ((int)reservoir.size() < RESERVOIR_SIZE) {
            reservoir.push_back({ phase, angle, c_b, c_r, s_b, s_r, e_c, e_s, mag });
            if ((int)reservoir.size() == RESERVOIR_SIZE) {
                std::make_heap(reservoir.begin(), reservoir.end(),
                    [](const SampleErr& a, const SampleErr& b) { return a.err_mag > b.err_mag; });
            }
        }
        else if (mag > reservoir.front().err_mag) {
            std::pop_heap(reservoir.begin(), reservoir.end(),
                [](const SampleErr& a, const SampleErr& b) { return a.err_mag > b.err_mag; });
            reservoir.back() = { phase, angle, c_b, c_r, s_b, s_r, e_c, e_s, mag };
            std::push_heap(reservoir.begin(), reservoir.end(),
                [](const SampleErr& a, const SampleErr& b) { return a.err_mag > b.err_mag; });
        }
    }

    const double rms_cos = std::sqrt(sum_sq_cos / (double)total);
    const double rms_sin = std::sqrt(sum_sq_sin / (double)total);

    // sort and output worst samples
    std::sort(reservoir.begin(), reservoir.end(),
        [](const SampleErr& a, const SampleErr& b) { return a.err_mag > b.err_mag; });
    const int top_n = std::min<int>(TOP_KEEP, (int)reservoir.size());

    {
        std::ostringstream oss;
        oss << "seam_scan_top_" << seam.label << ".csv";
        std::ofstream topCsv(oss.str().c_str(), std::ios::binary);
        topCsv << std::setprecision(17) << std::fixed;
        topCsv << "rank,phase,angle_rad,angle_deg,cos_btrig,cos_std,cos_abs_err,"
            "sin_btrig,sin_std,sin_abs_err,err_mag\n";
        for (int i = 0; i < top_n; ++i) {
            const auto& s = reservoir[i];
            const double deg = s.phase * 360.0;
            topCsv << (i + 1) << "," << s.phase << "," << s.angle << "," << deg << ","
                << s.cos_b << "," << s.cos_r << "," << s.err_c_abs << ","
                << s.sin_b << "," << s.sin_r << "," << s.err_s_abs << ","
                << s.err_mag << "\n";
        }
    }

    // pass 2: histogram
    std::vector<long long> hist(HIST_BINS, 0);
    auto bin_index = [&](double mag) -> int {
        if (global_max_mag <= 0.0) return 0;
        int b = (int)std::floor(mag / global_max_mag * (HIST_BINS - 1));
        if (b < 0) b = 0; if (b >= HIST_BINS) b = HIST_BINS - 1;
        return b;
        };

    for (long long k = -STEPS_PER_SIDE; k <= STEPS_PER_SIDE; ++k) {
        const double dk = static_cast<double>(k);
        const double phase = std::fma(dk, PHASE_STEP, center_phase);    // CHANGE: no drift
        const double angle = std::fma(phase, btrig::TAU, 0.0);          // CHANGE: fused scale

        double c_b, s_b;
        btrig::sincos(phase, USE_PRECISE, c_b, s_b);
        const double c_r = std::cos(angle);
        const double s_r = std::sin(angle);

        const double e_c = std::fabs(c_b - c_r);
        const double e_s = std::fabs(s_b - s_r);
        const double mag = std::sqrt(e_c * e_c + e_s * e_s);

        hist[bin_index(mag)]++;
    }

    auto t1 = Clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    std::cout << "RMS abs error (cos): " << rms_cos << "\n";
    std::cout << "RMS abs error (sin): " << rms_sin << "\n";
    std::cout << "Max abs error (cos): " << max_cos_err
        << " at phase " << max_cos_err_phase
        << " (angle " << (max_cos_err_phase * btrig::TAU)
        << " rad, " << (max_cos_err_phase * 360.0) << " deg)\n";
    std::cout << "Max abs error (sin): " << max_sin_err
        << " at phase " << max_sin_err_phase
        << " (angle " << (max_sin_err_phase * btrig::TAU)
        << " rad, " << (max_sin_err_phase * 360.0) << " deg)\n";
    std::cout << "Global max combined error magnitude: " << global_max_mag << "\n";
    std::cout << "Elapsed: " << ms << " ms\n";

    {
        std::ostringstream ossTop, ossHist;
        ossTop << "seam_scan_top_" << seam.label << ".csv";
        ossHist << "seam_scan_hist_" << seam.label << ".csv";
        std::cout << "Wrote top " << top_n << " to: " << ossTop.str() << "\n";
        std::cout << "Wrote histogram to: " << ossHist.str() << "\n";
    }
}


int main() {
    std::cout.setf(std::ios::fixed); std::cout << std::setprecision(17);

    // Print the list of tested seams up front for clarity
    std::cout << "Testing known difficult sin/cos seams:\n";
    for (int i = 0; i < NUM_SEAMS; ++i) {
        const double phase = SEAMS[i].phase;
        std::cout << "  - " << SEAMS[i].label
            << " | phase=" << phase
            << " | angle_rad=" << (phase * btrig::TAU)
            << " | angle_deg=" << (phase * 360.0)
            << "\n";
    }

    for (int i = 0; i < NUM_SEAMS; ++i) {
        scan_seam(SEAMS[i]);
    }

    // ---------------------------------
    // Throughput benchmark
    // ---------------------------------
    std::cout.setf(std::ios::fixed); std::cout << std::setprecision(6);
    std::cout << "Generating " << NUM_SAMPLES << " phases...\n";
    auto phases = make_phases(NUM_SAMPLES);

    // Warm-up (ICACHE)
    {
        double sink = 0.0;
        for (size_t i = 0; i < 10000 && i < phases.size(); ++i) {
            double c, s;
            btrig::sincos(phases[i], USE_PRECISE, c, s);
            sink += c - s * 0.123;
        }
        volatile double waste = sink; (void)waste;
    }

    // Benchmark: std::cos + std::sin (radian ref)
    double std_time_ms = 0.0, std_sink = 0.0;
    for (int rep = 0; rep < REPEATS; ++rep) {
        double start = now_ms();
        double sink = 0.0;
        for (size_t i = 0; i < phases.size(); ++i) {
            const double angle = phases[i] * btrig::TAU;
            const double c = std::cos(angle);
            const double s = std::sin(angle);
            sink += c * 0.5 + s * 0.25;
        }
        double end = now_ms();
        std_time_ms += (end - start);
        std_sink += sink;
    }

    // Benchmark: btrig::sincos (phase pair)
    double pair_time_ms = 0.0, pair_sink = 0.0;
    for (int rep = 0; rep < REPEATS; ++rep) {
        double start = now_ms();
        double sink = 0.0;
        for (size_t i = 0; i < phases.size(); ++i) {
            double c, s;
            btrig::sincos(phases[i], USE_PRECISE, c, s);
            sink += c * 0.5 + s * 0.25;
        }
        double end = now_ms();
        pair_time_ms += (end - start);
        pair_sink += sink;
    }

    // Optional: singletons (kept for visibility; you can remove)
    double cos_time_ms = 0.0, cos_sink = 0.0;
    for (int rep = 0; rep < REPEATS; ++rep) {
        double start = now_ms();
        double sink = 0.0;
        double c;
        for (size_t i = 0; i < phases.size(); ++i) {
            
            btrig::cos(phases[i], USE_PRECISE,c);
            sink += c;
        }
        double end = now_ms();
        cos_time_ms += (end - start);
        cos_sink += sink;
    }

    double sin_time_ms = 0.0, sin_sink = 0.0;
    for (int rep = 0; rep < REPEATS; ++rep) {
        double start = now_ms();
        double sink = 0.0;
        double s;
        for (size_t i = 0; i < phases.size(); ++i) {
            btrig::sin(phases[i], USE_PRECISE, s);
            sink += s;
        }
        double end = now_ms();
        sin_time_ms += (end - start);
        sin_sink += sink;
    }

    // Normalize by repeats
    std_time_ms /= REPEATS;
    pair_time_ms /= REPEATS;
    cos_time_ms /= REPEATS;
    sin_time_ms /= REPEATS;

    // Throughput (million evals/sec)
    const double n = static_cast<double>(NUM_SAMPLES);
    const double std_mpps = (n / std_time_ms) / 1e3;
    const double pair_mpps = (n / pair_time_ms) / 1e3;
    const double cos_mcps = (n / cos_time_ms) / 1e3;
    const double sin_mcps = (n / sin_time_ms) / 1e3;

    std::cout << "\n--- TRIG SPEED TEST (" << NUM_SAMPLES << " samples, phase-driven, "
        << (USE_PRECISE ? "precise" : "fast") << ") ---\n";

    std::cout << "std::sin/cos pair (calls)        : "
        << std_time_ms << " ms    (sink=" << std_sink << ")  ["
        << std_mpps << " M pairs/s]\n";

    std::cout << "btrig::sincos (phase pair)       : "
        << pair_time_ms << " ms    (sink=" << pair_sink << ")  ["
        << pair_mpps << " M pairs/s]\n";

    std::cout << "btrig::cos only (phase singleton): "
        << cos_time_ms << " ms    (sink=" << cos_sink << ")  ["
        << cos_mcps << " M calls/s]\n";

    std::cout << "btrig::sin only (phase singleton): "
        << sin_time_ms << " ms    (sink=" << sin_sink << ")  ["
        << sin_mcps << " M calls/s]\n";

    std::cout << std::setprecision(3);
    std::cout << "\nSpeedup vs std pair (pairs): "
        << (std_time_ms / pair_time_ms) << "x\n";

    return 0;
}
