#include <chrono>
#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include "btrig.h"

using namespace std;
using namespace btrig;

template<typename F>
double bench(const char* name, F func, int N = 5'000'000) {
    std::mt19937_64 rng(12345);
    std::uniform_real_distribution<double> dist(-TAU, TAU);

    volatile double sink = 0.0;
    for (int i = 0; i < 1000; ++i) sink += func(dist(rng), dist(rng));

    auto t0 = chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i) {
        double a = dist(rng), b = dist(rng);
        sink += func(a, b);
    }
    auto t1 = chrono::high_resolution_clock::now();

    double ms = chrono::duration<double, std::milli>(t1 - t0).count();
    cout << left << setw(28) << name
        << " : " << setw(10) << ms
        << " ms    (sink=" << sink << ")\n";
    return ms;
}

int main() {
    const int N = 5'000'000;

    cout << "\n--- TRIG SPEED TEST (" << N << " samples) ---\n\n";

    // sin/cos pair
    bench("std::sin/cos pair", [](double a, double) {
        return std::sin(a) + std::cos(a);
        }, N);
    bench("btrig::sincos fast", [](double a, double) {
        double c, s; btrig::sincos(a, false, c, s);
        return c + s;
        }, N);
    bench("btrig::sincos precise", [](double a, double) {
        double c, s; btrig::sincos(a, true, c, s);
        return c + s;
        }, N);

    cout << "\n";

    // singletons
    bench("std::sin only", [](double a, double) {
        return std::sin(a);
        }, N);
    bench("btrig::sin fast", [](double a, double) {
        double s; btrig::sin(a, false, s); return s;
        }, N);
    bench("btrig::sin precise", [](double a, double) {
        double s; btrig::sin(a, true, s); return s;
        }, N);

    cout << "\n";

    bench("std::cos only", [](double a, double) {
        return std::cos(a);
        }, N);
    bench("btrig::cos fast", [](double a, double) {
        double c; btrig::cos(a, false, c); return c;
        }, N);
    bench("btrig::cos precise", [](double a, double) {
        double c; btrig::cos(a, true, c); return c;
        }, N);

    cout << "\n";

    // tan
    bench("std::tan", [](double a, double) {
        return std::tan(a);
        }, N);
    bench("btrig::tan", [](double a, double) {
        return btrig::tan(a);
        }, N);

    cout << "\n";

    // atan2
    bench("std::atan2", [](double y, double x) {
        return std::atan2(y, x);
        }, N);
    bench("btrig atan2", [](double y, double x) {
        return btrig::atan2(y, x);
        }, N);

    cout << "\n";

    // atan
    bench("std::atan", [](double y, double) {
        return std::atan(y);
        }, N);
    bench("btrig atan", [](double y, double) {
        return btrig::atan(y);
        }, N);

    cout << "\n";

    // asin
    bench("std::asin", [](double z, double) {
        return std::asin(std::max(-1.0, std::min(1.0, z)));
        }, N);
    bench("btrig asin", [](double z, double) {
        return btrig::asin(z);
        }, N);

    cout << "\n";

    // acos
    bench("std::acos", [](double z, double) {
        return std::acos(std::max(-1.0, std::min(1.0, z)));
        }, N);
    bench("btrig acos", [](double z, double) {
        return btrig::acos(z);
        }, N);

    cout << "\n";

    // symbolic angle
    bench("btrig angle from xy", [](double x, double y) {
        return btrig::atan2(x, y);
        }, N);

    cout << "\n--- END TEST ---\n";
    return 0;
}