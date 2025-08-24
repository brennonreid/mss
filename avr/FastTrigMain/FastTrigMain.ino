#include "btrig.h"
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
static inline float u01(uint64_t &state) {
    return (float)((splitmix64(state) >> 41) * (1.0f / (1ULL << 23)));
}

// ---------------- Config ----------------
const uint32_t numSamples = 10000;
constexpr float TAU = 6.2831853071795864769f;
constexpr float ONE_OVER_TAU = 1.0f / TAU;
const uint64_t SEED = 0xCAFEBABEDEADBEEFULL;

// Volatile sink to prevent compiler optimizations
volatile float sink = 0.0f;

// Helper to print benchmark results
void print_result(const char* name, unsigned long time_us, float checksum) {
    Serial.print(name);
    Serial.print(" : ");
    Serial.print(time_us);
    Serial.print(" us    (sink=");
    Serial.print(checksum, 8);
    Serial.println(")");
}

void setup() {
    Serial.begin(115200);
    while (!Serial) {}
    Serial.println("btrig:: lookup tables self-initialize...");
    delay(2000);
    Serial.println("\n--- TRIG SPEED TEST ---\n");
}

void loop() {
    uint64_t seed = SEED;

    // Bench: sin/cos pair
    Serial.println("--- sin/cos pairs ---");
    unsigned long start_time, end_time;
    float c, s;

    // std::sin/cos pair
    seed = SEED;
    sink = 0.0f;
    start_time = micros();
    for (uint32_t i = 0; i < numSamples; ++i) {
        float angle = u01(seed) * TAU;
        c = cos(angle);
        s = sin(angle);
        sink += c + s;
    }
    end_time = micros();
    print_result("std::sin/cos pair", end_time - start_time, sink);

    // btrig::sincos fast
    seed = SEED;
    sink = 0.0f;
    start_time = micros();
    for (uint32_t i = 0; i < numSamples; ++i) {
        float angle = u01(seed) * TAU;
        btrig::sincos(angle, false, c, s);
        sink += c + s;
    }
    end_time = micros();
    print_result("btrig::sincos fast", end_time - start_time, sink);
    Serial.println();

    // Bench: singletons
    Serial.println("--- singletons ---");
    float val;

    // std::sin only
    seed = SEED;
    sink = 0.0f;
    start_time = micros();
    for (uint32_t i = 0; i < numSamples; ++i) {
        float angle = u01(seed) * TAU;
        val = sin(angle);
        sink += val;
    }
    end_time = micros();
    print_result("std::sin only", end_time - start_time, sink);

    // btrig::sin fast
    seed = SEED;
    sink = 0.0f;
    start_time = micros();
    for (uint32_t i = 0; i < numSamples; ++i) {
        float angle = u01(seed) * TAU;
        btrig::sin(angle, false, val);
        sink += val;
    }
    end_time = micros();
    print_result("btrig::sin fast", end_time - start_time, sink);
    Serial.println();

    // std::cos only
    seed = SEED;
    sink = 0.0f;
    start_time = micros();
    for (uint32_t i = 0; i < numSamples; ++i) {
        float angle = u01(seed) * TAU;
        val = cos(angle);
        sink += val;
    }
    end_time = micros();
    print_result("std::cos only", end_time - start_time, sink);

    // btrig::cos fast
    seed = SEED;
    sink = 0.0f;
    start_time = micros();
    for (uint32_t i = 0; i < numSamples; ++i) {
        float angle = u01(seed) * TAU;
        btrig::cos(angle, false, val);
        sink += val;
    }
    end_time = micros();
    print_result("btrig::cos fast", end_time - start_time, sink);
    Serial.println();
    
    // Bench: tan
    Serial.println("--- tangents ---");

    // std::tan
    seed = SEED;
    sink = 0.0f;
    start_time = micros();
    for (uint32_t i = 0; i < numSamples; ++i) {
        float angle = u01(seed) * TAU;
        val = tan(angle);
        sink += val;
    }
    end_time = micros();
    print_result("std::tan", end_time - start_time, sink);

    // btrig::tan
    seed = SEED;
    sink = 0.0f;
    start_time = micros();
    for (uint32_t i = 0; i < numSamples; ++i) {
        float angle = u01(seed) * TAU;
        val = btrig::tan(angle);
        sink += val;
    }
    end_time = micros();
    print_result("btrig::tan", end_time - start_time, sink);
    Serial.println();

    // Bench: atan2
    Serial.println("--- atan2 ---");
    float y, x;

    // std::atan2
    seed = SEED;
    sink = 0.0f;
    start_time = micros();
    for (uint32_t i = 0; i < numSamples; ++i) {
        y = u01(seed) * 2.0f - 1.0f;
        x = u01(seed) * 2.0f - 1.0f;
        val = atan2(y, x);
        sink += val;
    }
    end_time = micros();
    print_result("std::atan2", end_time - start_time, sink);

    // btrig::atan2
    seed = SEED;
    sink = 0.0f;
    start_time = micros();
    for (uint32_t i = 0; i < numSamples; ++i) {
        y = u01(seed) * 2.0f - 1.0f;
        x = u01(seed) * 2.0f - 1.0f;
        val = btrig::atan2(y, x, true);
        sink += val;
    }
    end_time = micros();
    print_result("btrig::atan2", end_time - start_time, sink);
    Serial.println();

    // Bench: atan
    Serial.println("--- atan ---");

    // std::atan
    seed = SEED;
    sink = 0.0f;
    start_time = micros();
    for (uint32_t i = 0; i < numSamples; ++i) {
        float angle = u01(seed) * TAU;
        val = atan(angle);
        sink += val;
    }
    end_time = micros();
    print_result("std::atan", end_time - start_time, sink);

    // btrig::atan
    seed = SEED;
    sink = 0.0f;
    start_time = micros();
    for (uint32_t i = 0; i < numSamples; ++i) {
        float angle = u01(seed) * TAU;
        val = btrig::atan(angle);
        sink += val;
    }
    end_time = micros();
    print_result("btrig::atan", end_time - start_time, sink);
    Serial.println();

    // Bench: asin
    Serial.println("--- asin ---");

    // std::asin
    seed = SEED;
    sink = 0.0f;
    start_time = micros();
    for (uint32_t i = 0; i < numSamples; ++i) {
        float z = u01(seed) * 2.0f - 1.0f;
        val = asin(z);
        sink += val;
    }
    end_time = micros();
    print_result("std::asin", end_time - start_time, sink);

    // btrig::asin
    seed = SEED;
    sink = 0.0f;
    start_time = micros();
    for (uint32_t i = 0; i < numSamples; ++i) {
        float z = u01(seed) * 2.0f - 1.0f;
        val = btrig::asin(z);
        sink += val;
    }
    end_time = micros();
    print_result("btrig::asin", end_time - start_time, sink);
    Serial.println();

    // Bench: acos
    Serial.println("--- acos ---");

    // std::acos
    seed = SEED;
    sink = 0.0f;
    start_time = micros();
    for (uint32_t i = 0; i < numSamples; ++i) {
        float z = u01(seed) * 2.0f - 1.0f;
        val = acos(z);
        sink += val;
    }
    end_time = micros();
    print_result("std::acos", end_time - start_time, sink);

    // btrig::acos
    seed = SEED;
    sink = 0.0f;
    start_time = micros();
    for (uint32_t i = 0; i < numSamples; ++i) {
        float z = u01(seed) * 2.0f - 1.0f;
        val = btrig::acos(z);
        sink += val;
    }
    end_time = micros();
    print_result("btrig::acos", end_time - start_time, sink);
    Serial.println();

    Serial.println("\n--- END TEST ---");
    delay(25000);
}