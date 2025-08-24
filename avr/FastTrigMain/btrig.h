#pragma once
#ifndef BTRIG_H
#define BTRIG_H

#include <math.h> // Required for fmodf

namespace btrig {

    // --------------------
    // Configuration
    // --------------------
    constexpr int TRIG_ANCHORS_BASE_POW = 4;
    constexpr int NUM_ANCHORS_QUADRANT = 1 << TRIG_ANCHORS_BASE_POW;
    constexpr int NUM_ANCHORS_FULL = NUM_ANCHORS_QUADRANT * 4;
    constexpr int ANCHOR_MASK = NUM_ANCHORS_FULL - 1;

    // --------------------
    // Core constants (changed to float)
    // --------------------
    constexpr float TAU = 6.283185307179586476925286766559f; // 2*pi
    constexpr float QUADTAU = TAU / 4.0f;
    constexpr float STEP = TAU / static_cast<float>(NUM_ANCHORS_FULL);
    constexpr float ONE_OVER_STEP = 1.0f / STEP;
    constexpr float INV_TAU = 1.0f / TAU;

    // Precomputed inverse factorials for interpolation
    constexpr float INV_2 = 1.0f / 2.0f;
    constexpr float INV_6 = 1.0f / 6.0f;
    constexpr float INV_24 = 1.0f / 24.0f;
    constexpr float INV_120 = 1.0f / 120.0f;
    constexpr float INV_720 = 1.0f / 720.0f;
    constexpr float INV_5040 = 1.0f / 5040.0f;
    constexpr float INV_40320 = 1.0f / 40320.0f;

    static float fullCircleAnchors[NUM_ANCHORS_FULL][3];

    // ---- Internal helpers ----
    static inline void applyQuadrantTransform(float inx, float iny, int quadrant,
        float& outx, float& outy) {
        switch (quadrant & 3) {
        case 0: outx = inx; outy = iny; break;
        case 1: outx = -iny; outy = inx; break;
        case 2: outx = -inx; outy = -iny; break;
        default: outx = iny; outy = -inx; break;
        }
    }
    
    static inline float wrapAngle(float angle) {
        float wrapped = fmodf(angle, TAU);
        if (wrapped < 0.0f) {
            wrapped += TAU;
        }
        return wrapped;
    }

    static inline float taylorCosFloat(float x, int depth) {
        float term = 1.0f;
        float sum = 1.0f;
        for (int i = 1; i <= depth; ++i) {
            term *= -x * x / ((2.0f * i - 1.0f) * (2.0f * i));
            sum += term;
        }
        return sum;
    }

    static inline float taylorSinFloat(float x, int depth) {
        float term = x;
        float sum = x;
        for (int i = 1; i <= depth; ++i) {
            term *= -x * x / ((2.0f * i) * (2.0f * i + 1.0f));
            sum += term;
        }
        return sum;
    }

    static inline void taylorSinCosFloat2(float angle, int depth,
        float& outSin, float& outCos) {
        angle = wrapAngle(angle);
        int quadrant = static_cast<int>(angle / QUADTAU);
        float x = angle - quadrant * QUADTAU;

        float sinVal = taylorSinFloat(x, depth);
        float cosVal = taylorCosFloat(x, depth);

        applyQuadrantTransform(cosVal, sinVal, quadrant, outCos, outSin);
    }

    static void initFullCircleAnchors() {
        for (int i = 0; i < NUM_ANCHORS_FULL; ++i) {
            float angle = i * STEP;

            int quadrant = static_cast<int>(angle / QUADTAU);
            float x = angle - quadrant * QUADTAU;

            float cosVal = taylorCosFloat(x, 10);
            float sinVal = taylorSinFloat(x, 10);

            float cosFinal, sinFinal;
            applyQuadrantTransform(cosVal, sinVal, quadrant, cosFinal, sinFinal);

            fullCircleAnchors[i][0] = angle;
            fullCircleAnchors[i][1] = cosFinal;
            fullCircleAnchors[i][2] = sinFinal;
        }
    }

    static struct _InitFullCircleAnchors {
        _InitFullCircleAnchors() { initFullCircleAnchors(); }
    } _initFullCircleAnchorsInstance;

    static inline void sin(float angle, bool precise, float& outSin) {
        const float t = wrapAngle(angle) * ONE_OVER_STEP;
        const int idx = static_cast<int>(t) & ANCHOR_MASK;
        const float d = (t - static_cast<int>(t)) * STEP;
        const float d2 = d * d;
        const float* a = fullCircleAnchors[idx];
        float dx, dy;
        if (precise) {
    dx = 1.0f - d2 * (INV_2 - d2 * INV_24);
    dy = d * (1.0f - d2 * INV_6);
} else {
            dx = 1.0f - d2 * INV_2;
            dy = d;
        }
        outSin = dx * a[2] + dy * a[1];
    }

    static inline void cos(float angle, bool precise, float& outCos) {
        const float t = wrapAngle(angle) * ONE_OVER_STEP;
        const int idx = static_cast<int>(t) & ANCHOR_MASK;
        const float d = (t - static_cast<int>(t)) * STEP;
        const float d2 = d * d;
        const float* a = fullCircleAnchors[idx];
        float dx, dy;
      if (precise) {
    dx = 1.0f - d2 * (INV_2 - d2 * INV_24);
    dy = d * (1.0f - d2 * INV_6);
} else {
            dx = 1.0f - d2 * INV_2;
            dy = d;
        }
        outCos = dx * a[1] - dy * a[2];
    }

    static inline void sincos(float angle, bool precise,
        float& cosOut, float& sinOut) {
        const float t = wrapAngle(angle) * ONE_OVER_STEP;
        const int idx = static_cast<int>(t) & ANCHOR_MASK;
        const float d = (t - static_cast<int>(t)) * STEP;
        const float d2 = d * d;
        const float* a = fullCircleAnchors[idx];
        float dx, dy;
       if (precise) {
    dx = 1.0f - d2 * (INV_2 - d2 * INV_24);
    dy = d * (1.0f - d2 * INV_6);
} else {
            dx = 1.0f - d2 * INV_2;
            dy = d;
        }
        cosOut = dx * a[1] - dy * a[2];
        sinOut = dx * a[2] + dy * a[1];
    }

    static inline float atan2(float y, float x, bool precise) {
    if (x == 0.0f && y == 0.0f) return 0.0f;

    const float ax = (y < 0.0f ? -y : y);
    const float ay = (x < 0.0f ? -x : x);

    const bool swap = (ax > ay);
    const float num = swap ? ay : ax;
    const float den = swap ? ax : ay;
    const float z = num / den;
    const float z2 = z * z;

    float a;
    if (precise) {
        a = z * (0.9998660f + z2 * (-0.3302995f + z2 * (0.1801410f + z2 * (-0.0851330f + 0.0208351f * z2))));
    } else {
        a = z * (0.9998660f + z2 * (-0.3302995f + z2 * (0.1801410f)));
    }

    const float angle = swap ? (1.5707963267948966f - a) : a;

    if (x >= 0.0f) {
        return (y >= 0.0f) ? angle : -angle;
    } else {
        return (y >= 0.0f) ? (3.141592653589793f - angle) : (angle - 3.141592653589793f);
    }
}
    static inline float fast_sqrt(float x) {
        if (x <= 0.0f) return 0.0f;
        union { float d; unsigned long u; } v;
        v.d = x;
        v.u = 0x5f375a86 - (v.u >> 1);
        float y = v.d;
        const float xhalf = 0.5f * x;
        y = y * (1.5f - xhalf * y * y);
        y = y * (1.5f - xhalf * y * y);
        return x * y;
    }

    static inline float asin(float z) {
        const float HALFPI = TAU * 0.25f;
        if (z >= 1.0f) return HALFPI;
        if (z <= -1.0f) return -HALFPI;

        const float t = 1.0f - z * z;
        const float c = (t > 0.0f) ? fast_sqrt(t) : 0.0f;
        // FIX: Added the precise flag to the atan2 call
        return atan2(z, c, /*precise=*/false);
    }

    static inline float acos(float x) {
        if (x >= 1.0f) return 0.0f;
        if (x <= -1.0f) return TAU * 0.5f;

        float t = 1.0f - x * x;
        float s = (t > 0.0f) ? fast_sqrt(t) : 0.0f;
        // FIX: Added the precise flag to the atan2 call
        return atan2(s, x, /*precise=*/false);
    }

static inline float tan(float angle) {
    // constants
        const float PI_2   = 1.57079632679489661923f;
    const float PI_4   = 0.78539816339744830962f;
    const float INV_PI = 0.31830988618379067154f;

    // reduce to [-PI/2, PI/2] using k = round(angle / PI)
    float kf = angle * INV_PI;
    int   k  = (int)(kf + (kf >= 0.0f ? 0.5f : -0.5f));
    float y  = angle - (float)k * PI;

    float ay = y >= 0.0f ? y : -y;

    // tan polynomial on [-PI/4, PI/4]:
    // tan t ≈ t + t^3*(1/3 + x*(2/15 + x*(17/315))) with x = t^2
    auto tan_poly = [](float t) {
        const float C1 = 0.3333333333333333f;   // 1/3
        const float C2 = 0.1333333333333333f;   // 2/15
        const float C3 = 0.05396825396825397f;  // 17/315
        float x = t * t;
        float p = (C3 * x + C2);
        p = (p * x + C1);
        return t + t * x * p;
    };

    if (ay <= PI_4) {
        return tan_poly(y);
    } else {
        // outer octant: tan(y) = 1 / tan(PI/2 - |y|)
        float z = PI_2 - ay;     // in [0, PI/4]
        float t = tan_poly(z);
        float r = 1.0f / t;
        return y >= 0.0f ? r : -r;
    }
}



    static inline float atan(float y) {
        // FIX: Added the precise flag to the atan2 call
        return atan2(y, 1.0f, /*precise=*/false);
    }

} // namespace btrig

#endif // BTRIG_H