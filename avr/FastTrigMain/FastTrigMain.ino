#include "TaylorTrig.h"
#include <stdint.h>
#include <math.h>

// ---------------- PRNG: SplitMix64 ----------------
static inline uint64_t splitmix64(uint64_t &x) {
    uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

// Uniform [0,1)
static inline double u01(uint64_t &state) {
    // Take top 53 bits -> double in [0,1)
    return ( (splitmix64(state) >> 11) * (1.0 / (1ULL << 53)) );
}

// ---------------- Config ----------------
const uint32_t numSamples = 10000;           // set as high as you like
constexpr double TAU = 6.2831853071795864769;

// Fixed seeds so both passes see identical angles
const uint64_t SEED_BASE  = 0xCAFEBABEDEADBEEFULL;
const uint64_t SEED_PRINT = 0x123456789ABCDEF0ULL;  // for the 10 sample prints

void setup() {
    Serial.begin(115200);
    while (!Serial) {}

    Serial.println("Initializing TaylorTrig lookup tables...");
    TaylorTrig::setupTaylorTrig();
    Serial.println("Initialization complete.");
}

void loop() {
    bool precise = false;

    float  x_approx, y_approx;   // CHANGED: was double
    double x_native, y_native;

    volatile double customChecksum = 0.0;
    volatile double nativeChecksum = 0.0;

    // ========== SPEED TEST ==========
    // Custom TaylorTrig Benchmark
    uint64_t seed = SEED_BASE; // reset seed so both runs see same angles
    unsigned long startCustom = micros();
    for (uint32_t i = 0; i < numSamples; i++) {
        double angle = u01(seed) * TAU;
        TaylorTrig::getUnitVectorFromAngle((float)angle, x_approx, y_approx, precise); // cast angle to float
        customChecksum += (double)x_approx + (double)y_approx; // CHANGED: cast to double for checksum
    }
    unsigned long endCustom = micros();
    unsigned long timeCustom = endCustom - startCustom;

    // Native Math.h Benchmark
    seed = SEED_BASE; // replay same angles
    unsigned long startNative = micros();
    for (uint32_t i = 0; i < numSamples; i++) {
        double angle = u01(seed) * TAU;
        x_native = cos(angle);
        y_native = sin(angle);
        nativeChecksum += x_native + y_native;
    }
    unsigned long endNative = micros();
    unsigned long timeNative = endNative - startNative;

    // Print Results (Speed)
    Serial.println("--- SPEED TEST ---");
    Serial.print("Total Calculations: ");
    Serial.println((unsigned long)numSamples);
    Serial.print("Custom TaylorTrig time: ");
    Serial.print(timeCustom);
    Serial.println(" microseconds");
    Serial.print("Native Math.h time:     ");
    Serial.print(timeNative);
    Serial.println(" microseconds");
    Serial.print("Custom Checksum: ");
    Serial.println(customChecksum, 12);
    Serial.print("Approx Checksum: ");
    Serial.println(nativeChecksum, 12);
    Serial.println();

    // ========== ACCURACY TEST ==========
    // Reset checksums so this section reflects only the accuracy pass
    customChecksum = 0.0;
    nativeChecksum = 0.0;

    seed = SEED_BASE; // replay again for accuracy totals
    for (uint32_t i = 0; i < numSamples; i++) {
        double angle = u01(seed) * TAU;

        TaylorTrig::getUnitVectorFromAngle((float)angle, x_approx, y_approx, precise); // cast angle
        x_native = cos(angle);
        y_native = sin(angle);

        customChecksum += (double)x_approx + (double)y_approx; // cast to double
        nativeChecksum += x_native + y_native;
    }

    Serial.println("--- ACCURACY TEST ---");
    Serial.print("Total Calculations: ");
    Serial.println((unsigned long)numSamples);
    Serial.print("Custom Checksum: ");
    Serial.println(customChecksum, 12);
    Serial.print("Approx Checksum: ");
    Serial.println(nativeChecksum, 12);
    Serial.println();

    // Print 10 sample comparisons (deterministic, separate seed)
    uint64_t printSeed = SEED_PRINT;
    for (int i = 0; i < 10; i++) {
        double angle = u01(printSeed) * TAU;

        TaylorTrig::getUnitVectorFromAngle((float)angle, x_approx, y_approx, precise); // cast angle
        x_native = cos(angle);
        y_native = sin(angle);

        Serial.print("Angle: ");
        Serial.print(angle, 6);
        Serial.print(" | Custom (x, y): (");
        Serial.print((double)x_approx, 12); // CHANGED: cast for consistent print
        Serial.print(", ");
        Serial.print((double)y_approx, 12); // CHANGED: cast for consistent print
        Serial.print(") | Approx (x, y): (");
        Serial.print(x_native, 12);
        Serial.print(", ");
        Serial.print(y_native, 12);
        Serial.println(")");
    }
    Serial.println();

    delay(25000);
}
