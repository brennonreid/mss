#pragma once
#include <stdint.h> // For uint32_t (used in NUM_ANCHORS)

// Fixed sizes (no overrides)
#define NUM_ANCHORS_QUADRANT 32
#define NUM_ANCHORS_FULL (4 * NUM_ANCHORS_QUADRANT)

namespace TaylorTrig {

// Struct to hold angle, cosine, and sine values for lookup table anchors.
// These must be float to match the calculations in TaylorTrig.cpp
struct Anchor {
    float theta; // Angle in radians
    float cos;   // Cosine value
    float sin;   // Sine value
};

// Arrays are defined in the .cpp with these exact sizes
extern Anchor canonicalAnchors[NUM_ANCHORS_QUADRANT];
extern Anchor fullCircleAnchors[NUM_ANCHORS_FULL];

// API (names and types must match definitions in .cpp)
// These functions operate on float values
float taylorCosFloat(float x, int depth);
float taylorSinFloat(float x, int depth);

// Function to initialize the TaylorTrig lookup tables.
// This must be called once, typically in setup().
void setupTaylorTrig();

// Hot path (float math) - This signature now matches your TaylorTrig.cpp implementation
// angleRadians: The input angle in radians (float)
// x: Reference to a float to store the calculated cosine (x-component)
// y: Reference to a float to store the calculated sine (y-component)
// precise: A boolean flag to switch between approximation precisions
void getUnitVectorFromAngle(float angleRadians, float& x, float& y, bool precise);

} // namespace TaylorTrig
