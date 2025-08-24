#pragma once
#include <vector>
#include <cmath>
#include "btrig.h"  // your trig source

namespace brot {

    // --- Type definitions ---
    using Vec3 = std::vector<double>;

    struct EulerAngles { double yaw, pitch, roll; };
    struct AxisAngle { double kx, ky, kz, theta; };

    struct SinCosResult { double cos, sin; };

    // --- Compact inline cos/sin oracle (unit-safe) ---
    inline SinCosResult cosSinFromAngle(double theta, bool precise = false) {
        // Normalize with wrapTau so we never stress large magnitudes
        theta = btrig::wrapTau(theta);
        double c, s;
        btrig::sincos(theta, precise, c, s);
        return { c, s };
    }

    // --- Direct axis–angle rotation (Rodrigues) ---
    inline Vec3 rotateAxisAngle(double kx, double ky, double kz,
        double theta, double x, double y, double z,
        bool precise = false)
    {
        theta = btrig::wrapTau(theta);
        double c, s; btrig::sincos(theta, precise, c, s);

        const double cx = ky * z - kz * y;
        const double cy = kz * x - kx * z;
        const double cz = kx * y - ky * x;

        const double kd = kx * x + ky * y + kz * z;
        const double oneMinusC = 1.0 - c;

        return {
            x * c + cx * s + kx * kd * oneMinusC,
            y * c + cy * s + ky * kd * oneMinusC,
            z * c + cz * s + kz * kd * oneMinusC
        };
    }

    // --- Axis helpers (unit axes only) ---
    inline Vec3 rotateX(double theta, double x, double y, double z) { return rotateAxisAngle(1, 0, 0, theta, x, y, z); }
    inline Vec3 rotateY(double theta, double x, double y, double z) { return rotateAxisAngle(0, 1, 0, theta, x, y, z); }
    inline Vec3 rotateZ(double theta, double x, double y, double z) { return rotateAxisAngle(0, 0, 1, theta, x, y, z); }

    // --- Chained rotate (Z → X → Y) ---
    inline Vec3 rotateZXY(double a, double b, double c, double x, double y, double z) {
        Vec3 v = rotateZ(a, x, y, z);
        v = rotateX(b, v[0], v[1], v[2]);
        v = rotateY(c, v[0], v[1], v[2]);
        return v;
    }

    // --- Euler readout (ZYX) from axis–angle ---
    inline EulerAngles eulerZYXFromAxisAngle(double kx, double ky, double kz, double theta) {
        theta = btrig::wrapTau(theta);
        const auto sc = cosSinFromAngle(theta);
        const double c = sc.cos, s = sc.sin, C = 1.0 - c;

        const double r00 = c + kx * kx * C;
        const double r01 = kx * ky * C - kz * s;
        const double r02 = kx * kz * C + ky * s;

        const double r10 = ky * kx * C + kz * s;
        const double r11 = c + ky * ky * C;
        const double r12 = ky * kz * C - kx * s;

        const double r20 = kz * kx * C - ky * s;
        const double r21 = kz * ky * C + kx * s;
        const double r22 = c + kz * kz * C;

        const double pitch = btrig::asin(-r20);
        const double cp = std::hypot(r00, r10);

        double yaw, roll;
        if (cp < 1e-6) {
            yaw = btrig::atan2(-r01, r11);
            roll = 0.0;
        }
        else {
            yaw = btrig::atan2(r10, r00);
            roll = btrig::atan2(r21, r22);
        }
        return { yaw, pitch, roll };
    }

    // --- Compose axis–angle → axis–angle (no matrix/quat exposed) ---
    inline AxisAngle composeAxisAngle(double ax, double ay, double az, double aTheta,
        double bx, double by, double bz, double bTheta)
    {
        // Normalize inputs to [0, TAU) — no wrap-π anywhere
        aTheta = btrig::wrapTau(aTheta);
        bTheta = btrig::wrapTau(bTheta);

        // Half-angles directly (avoids sqrt( (1±cos)/2 ) and sign logic)
        double cah, sah; btrig::sincos(0.5 * aTheta, /*precise=*/false, cah, sah);
        double cbh, sbh; btrig::sincos(0.5 * bTheta, /*precise=*/false, cbh, sbh);

        // “Quaternion-like” composition in minimal form
        const double vax = ax * sah, vay = ay * sah, vaz = az * sah;
        const double vbx = bx * sbh, vby = by * sbh, vbz = bz * sbh;

        const double w = cbh * cah - (vbx * vax + vby * vay + vbz * vaz);
        const double vx = cbh * vax + cah * vbx + (vay * vbz - vaz * vby);
        const double vy = cbh * vay + cah * vby + (vaz * vbx - vax * vbz);
        const double vz = cbh * vaz + cah * vbz + (vax * vby - vay * vbx);

        double sn = std::hypot(std::hypot(vx, vy), vz);
        if (sn < 1e-18) return { 1.0, 0.0, 0.0, 0.0 };   // identity

        // Angle ∈ (0, 2π); normalize to [0, TAU) explicitly (no wrap-π)
        double theta = 2.0 * btrig::atan2(sn, w);
        theta = btrig::wrapTau(theta);

        const double inv = 1.0 / sn;
        return { vx * inv, vy * inv, vz * inv, theta };
    }

    // --- Rotate via rotation vector ---
    inline Vec3 rotateVectorFromRotVec(double rx, double ry, double rz,
        double x, double y, double z)
    {
        const double theta = std::hypot(rx, std::hypot(ry, rz));
        if (theta == 0.0) return { x, y, z };
        const double inv = 1.0 / theta;
        return rotateAxisAngle(rx * inv, ry * inv, rz * inv, theta, x, y, z);
    }

    inline EulerAngles eulerZYXFromRotVec(double rx, double ry, double rz) {
        const double theta = std::hypot(rx, std::hypot(ry, rz));
        if (theta == 0.0) return { 0.0, 0.0, 0.0 };
        const double inv = 1.0 / theta;
        return eulerZYXFromAxisAngle(rx * inv, ry * inv, rz * inv, theta);
    }

} // namespace brot
