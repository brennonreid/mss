#pragma once
#ifndef TAYLORTRIG_H
#define TAYLORTRIG_H

#include <vector>
#include <cmath>

struct Anchor {
    double theta;
    double cos;
    double sin;
};

// ---- Constants (exact JS equivalents) ----
extern const double TAU;
extern const double QUADTAU;

extern const int    NUM_ANCHORS_QUADRANT;
extern const int    NUM_ANCHORS_FULL;

extern const double DEG_STEP_QUADRANT;
extern const double DEG_STEP_FULL;
extern const double STEP;
extern const double ONE_OVER_STEP;

// Precomputed inverse factorials
extern const double INV_2;
extern const double INV_6;
extern const double INV_24;
extern const double INV_120;
extern const double INV_720;
extern const double INV_5040;
extern const double INV_40320;
extern const double INV_TAU;

// ---- Anchor storage (made external so header-only funcs can use them) ----
extern std::vector<Anchor> canonicalAnchors;
extern std::vector<Anchor> fullCircleAnchors;

// ---- Misc public helpers (unchanged) ----
double wrapTau(double a);
double taylorCosFloat(double x, int depth = 2);
double taylorSinFloat(double x, int depth = 2);

double angle_from_xy_symbolic(double x, double y, int M = 64, int refine = 18);
double atan2_sym(double y, double x);
double asin_sym(double z);

// ======================================================================
// Header-only hot path: getUnitVectorFromAngle2 / customSin / customCos
// ======================================================================

// Vector result via out-params
static __forceinline void getUnitVectorFromAngle2(double angle, bool precise,
    double& cosOut, double& sinOut)
{
    if (angle >= TAU || angle < 0.0) { angle -= std::floor(angle * INV_TAU) * TAU; }

    /*
    if (angle >= TAU || angle < 0.0) {
        angle = std::fmod(angle, TAU);
        if (angle < 0.0) angle += TAU;
    }
    */

    int idx = (int)(angle * ONE_OVER_STEP);
    const Anchor& a = fullCircleAnchors[(size_t)idx];
    const double d = angle - a.theta;

    const double d2 = d * d;

    const double dx = precise
        ? 1.0 - d2 * (INV_2 - d2 * (INV_24 - d2 * (INV_720 - d2 * INV_40320)))
        : 1.0 - d * d * INV_2;

    const double dy = precise
        ? d * (1.0 - d2 * (INV_6 - d2 * (INV_120 - d2 * INV_5040)))
        : d;

    cosOut = dx * a.cos - dy * a.sin;
    sinOut = dx * a.sin + dy * a.cos;
}

// Cos-only
static __forceinline void customCos(double angle, bool precise, double& cOut)
{
    if (angle >= TAU || angle < 0.0) { angle -= std::floor(angle * INV_TAU) * TAU; }

    int idx = (int)(angle * ONE_OVER_STEP);
    const Anchor& a = fullCircleAnchors[(size_t)idx];
    const double d = angle - a.theta;

    const double d2 = d * d;

    const double dx = precise
        ? 1.0 - d2 * (INV_2 - d2 * (INV_24 - d2 * (INV_720 - d2 * INV_40320)))
        : 1.0 - d2 * INV_2;  // reuse d2

    const double dy = precise
        ? d * (1.0 - d2 * (INV_6 - d2 * (INV_120 - d2 * INV_5040)))
        : d;

    // cosine rotation (differs from customSin)
    cOut = dx * a.cos - dy * a.sin;
}



static __forceinline void customSin(double angle, bool precise, double& sOut) {
    if (angle >= TAU || angle < 0.0) { angle -= std::floor(angle * INV_TAU) * TAU; }

    int idx = (int)(angle * ONE_OVER_STEP);
    const Anchor& a = fullCircleAnchors[(size_t)idx];
    const double d = angle - a.theta;

    const double d2 = d * d;

    const double dx = precise
        ? 1.0 - d2 * (INV_2 - d2 * (INV_24 - d2 * (INV_720 - d2 * INV_40320)))
        : 1.0 - d * d * INV_2;

    const double dy = precise
        ? d * (1.0 - d2 * (INV_6 - d2 * (INV_120 - d2 * INV_5040)))
        : d;

    sOut = dx * a.sin + dy * a.cos;
}


#endif // TAYLORTRIG_H
