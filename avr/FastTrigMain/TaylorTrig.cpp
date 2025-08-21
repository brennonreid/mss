#include "TaylorTrig.h"
#include <math.h>

namespace TaylorTrig {

// Constants
static constexpr float TAU      = 6.283185307179586f;
static constexpr float QUAD_TAU = TAU / 4.0f;

static constexpr float DEG_STEP_QUADRANT = QUAD_TAU / float(NUM_ANCHORS_QUADRANT);
static constexpr float DEG_STEP_FULL     = TAU / float(NUM_ANCHORS_FULL);
static constexpr float INV_DEG_STEP_FULL = 1.0f / DEG_STEP_FULL;

// Storage
Anchor canonicalAnchors[NUM_ANCHORS_QUADRANT];
Anchor fullCircleAnchors[NUM_ANCHORS_FULL];

// High-depth Taylor (no wrapping; caller supplies x)
float taylorCosFloat(float x, int depth) {
    float x2   = x * x;
    float term = 1.0f;  // k = 0
    float y    = term;
    for (int k = 1; k <= depth; ++k) {
        term *= -x2 / ((2.0f * k - 1.0f) * (2.0f * k));  // (-1)^k x^(2k)/(2k)!
        y    += term;
    }
    return y;
}

float taylorSinFloat(float x, int depth) {
    float x2   = x * x;
    float term = x;     // k = 0
    float y    = term;
    for (int k = 1; k <= depth; ++k) {
        term *= -x2 / ((2.0f * k) * (2.0f * k + 1.0f));  // (-1)^k x^(2k+1)/(2k+1)!
        y    += term;
    }
    return y;
}

// Build LUTs
void setupTaylorTrig() {
    // Canonical quadrant
    for (int i = 0; i < NUM_ANCHORS_QUADRANT; ++i) {
        const float angle = float(i) * DEG_STEP_QUADRANT;
        canonicalAnchors[i].theta = angle;
        canonicalAnchors[i].cos   = taylorCosFloat(angle, 55);
        canonicalAnchors[i].sin   = taylorSinFloat(angle, 55);
    }

    // Full circle via symmetry
    for (int i = 0; i < NUM_ANCHORS_FULL; ++i) {
        const float angle    = float(i) * DEG_STEP_FULL;
        const int quadrant   = int(floorf(angle / QUAD_TAU)) & 3;
        const float theta    = angle - float(quadrant) * QUAD_TAU;

        long idx = lroundf(theta / DEG_STEP_QUADRANT);
        if (idx >= NUM_ANCHORS_QUADRANT) idx = NUM_ANCHORS_QUADRANT - 1;

        const Anchor& c = canonicalAnchors[idx];
        Anchor a;
        switch (quadrant) {
            case 0: a = { angle,  c.cos,  c.sin }; break;
            case 1: a = { angle, -c.sin,  c.cos }; break;
            case 2: a = { angle, -c.cos, -c.sin }; break;
            default:a = { angle,  c.sin, -c.cos }; break;
        }
        fullCircleAnchors[i] = a;
    }
}

static constexpr float R1_2   = 1.0f / 2.0f;
static constexpr float R1_6   = 1.0f / 6.0f;
static constexpr float R1_24  = 1.0f / 24.0f;
static constexpr float R1_120 = 1.0f / 120.0f;

// Hot path (float)
void getUnitVectorFromAngle(float angleRadians, float& x, float& y, bool precise) {
    // Taylor coeffs
    
    // Normalize only if needed (fast path if already [0,TAU))
    float relativeAngle = angleRadians;

    // Fast index (truncate via multiply by inverse step)
    int zoneIndex = int(relativeAngle * INV_DEG_STEP_FULL) | 0;
    //if (zoneIndex >= NUM_ANCHORS_FULL) zoneIndex -= NUM_ANCHORS_FULL;

    const Anchor& anchor = fullCircleAnchors[zoneIndex];

    // Delta from LUT theta
    const float delta = relativeAngle - anchor.theta;
    const float d2    = delta * delta;

    float cosRel, sinRel;
    if (!precise) {
        // cos ≈ 1 - d^2/2 ;  sin ≈ d - d^3/6
        cosRel = 1.0f - d2 * R1_2;
        sinRel = delta * (1.0f - d2 * R1_6);
    } else {
        // Short Horner
        cosRel = 1.0f - d2 * (R1_2 - d2 * R1_24);
        sinRel = delta * (1.0f - d2 * (R1_6 - d2 * R1_120));
    }

    x = anchor.cos * cosRel - anchor.sin * sinRel;
    y = anchor.cos * sinRel + anchor.sin * cosRel;
}

} // namespace TaylorTrig
