#include "TaylorTrig.h"
#include <cmath>
#include <vector>

// ---- Constants ----
const double TAU = 6.283185307179586476925286766559;
const double QUADTAU = TAU / 4.0;

const int    NUM_ANCHORS_QUADRANT = 64;
const int    NUM_ANCHORS_FULL = 4 * NUM_ANCHORS_QUADRANT;

const double DEG_STEP_QUADRANT = QUADTAU / NUM_ANCHORS_QUADRANT;
const double DEG_STEP_FULL = TAU / NUM_ANCHORS_FULL;
const double STEP = TAU / NUM_ANCHORS_FULL;
const double ONE_OVER_STEP = 1.0 / STEP;

// Precomputed inverse factorials
const double INV_2 = 1.0 / 2.0;
const double INV_6 = 1.0 / 6.0;
const double INV_24 = 1.0 / 24.0;
const double INV_120 = 1.0 / 120.0;
const double INV_720 = 1.0 / 720.0;
const double INV_5040 = 1.0 / 5040.0;
const double INV_40320 = 1.0 / 40320.0;
const double INV_TAU = 1.0 / TAU;

// ---- Anchor tables ----
// ---- Anchor tables (quadrant via Taylor@55, full via symmetry) ----
std::vector<Anchor> canonicalAnchors = []() {
    std::vector<Anchor> anchors;
    anchors.reserve(NUM_ANCHORS_QUADRANT + 1);
    for (int i = 0; i <= NUM_ANCHORS_QUADRANT; ++i) {
        double theta = i * DEG_STEP_QUADRANT;
        // Use your Taylor functions at depth 55 for the LUT
        double c = taylorCosFloat(theta, 55);
        double s = taylorSinFloat(theta, 55);
        anchors.push_back({ theta, c, s });
    }
    return anchors;
    }();

std::vector<Anchor> fullCircleAnchors = []() {
    std::vector<Anchor> anchors;
    anchors.reserve(NUM_ANCHORS_FULL + 1);
    // Build by symmetry from the quadrant table
    for (int i = 0; i <= NUM_ANCHORS_FULL; ++i) {
        int q = i / NUM_ANCHORS_QUADRANT;      // 0..3 (quadrant)
        int t = i % NUM_ANCHORS_QUADRANT;      // 0..(NQ-1) position inside quadrant
        const Anchor& b = canonicalAnchors[(size_t)t];

        double c, s;
        switch (q & 3) {
        case 0: // [0, 90)
            c = b.cos;   s = b.sin;   break;
        case 1: // [90, 180)
            c = -b.sin;  s = b.cos;   break;
        case 2: // [180, 270)
            c = -b.cos;  s = -b.sin;  break;
        default: // [270, 360)
            c = b.sin;  s = -b.cos;  break;
        }

        double theta = i * STEP;
        anchors.push_back({ theta, c, s });
    }
    return anchors;
    }();

// ---- Utilities ----
double wrapTau(double a) {
    a = std::fmod(a, TAU);
    return (a < 0) ? a + TAU : a;
}

double taylorCosFloat(double x, int depth) {
    // cos(x) = sum_{k=0..depth-1} (-1)^k x^(2k) / (2k)!
    if (depth <= 1) return 1.0; // 1 term: 1
    const double x2 = x * x;
    double term = 1.0;  // k=0
    double sum = term;
    for (int k = 1; k < depth; ++k) {
        // term_k = term_{k-1} * ( -x^2 / ((2k-1)*(2k)) )
        term *= -x2 / ((double)(2 * k - 1) * (double)(2 * k));
        sum += term;
    }
    return sum;
}

double taylorSinFloat(double x, int depth) {
    // sin(x) = sum_{k=0..depth-1} (-1)^k x^(2k+1) / (2k+1)!
    if (depth <= 0) return 0.0; // 0 terms
    const double x2 = x * x;
    double term = x;   // k=0
    double sum = term;
    for (int k = 1; k < depth; ++k) {
        // term_k = term_{k-1} * ( -x^2 / ((2k)*(2k+1)) )
        term *= -x2 / ((double)(2 * k) * (double)(2 * k + 1));
        sum += term;
    }
    return sum;
}


double angle_from_xy_symbolic(double x, double y, int M, int refine) {
    double angle = 0.0;
    for (int i = 0; i < refine; ++i) {
        double minErr = 1e300;
        int best = 0;

        const double step = TAU / static_cast<double>(M);
        double theta = 0.0;
        for (int j = 0; j < M; ++j, theta += step) {
            double c, s;
            getUnitVectorFromAngle2(theta, /*precise=*/true, c, s);
            const double dx = c - x;
            const double dy = s - y;
            const double err = dx * dx + dy * dy;
            if (err < minErr) { minErr = err; best = j; }
        }
        angle = best * step;
        M <<= 1; // refine by doubling samples
    }
    return angle;
}


double atan2_sym(double y, double x) {
    if (x == 0.0 && y == 0.0) return 0.0;
    double bestErr = 1e300;
    int best = 0;
    for (int j = 0; j < NUM_ANCHORS_FULL; j++) {
        double dx = fullCircleAnchors[j].cos - x;
        double dy = fullCircleAnchors[j].sin - y;
        double err = dx * dx + dy * dy;
        if (err < bestErr) { bestErr = err; best = j; }
    }
    return fullCircleAnchors[best].theta;
}

double asin_sym(double z) {
    if (z > 1.0) z = 1.0;
    if (z < -1.0) z = -1.0;
    double bestErr = 1e300;
    int best = 0;
    for (int j = 0; j < NUM_ANCHORS_FULL; j++) {
        double dy = fullCircleAnchors[j].sin - z;
        double err = dy * dy;
        if (err < bestErr) { bestErr = err; best = j; }
    }
    return fullCircleAnchors[best].theta;
}
