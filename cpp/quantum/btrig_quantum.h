#pragma once
#ifndef BTRIG_QUANTUM_H
#define BTRIG_QUANTUM_H

// ----------------------------------------------------------------------------
// btrig_quantum.h — enhanced
// ----------------------------------------------------------------------------
// A lightweight, header-only companion for btrig focused on *turns*-normalized
// phases and practical quantum program emission (OpenQASM 3).
//
// What’s new in this enhanced version
//  • Exact dyadic emission: emit RZ/RX/RY as 2*pi*k/2^n (no rounding).
//  • Angle modes: choose Float64 or Dyadic per call.
//  • QFT-friendly phase schedules.
//  • (Optional) DDS streamer that batches sincos via btrig_simd for pulse-ish use.
//  • Fix: PhaseQ::half() now truly adds 1/2 turn, not 1 LSB.
//
// Build: C++17. This header stays portability-light on purpose.
// ----------------------------------------------------------------------------

#include <cstdint>
#include <cmath>
#include <string>
#include <vector>
#include <sstream>
#include <type_traits>

// If you want the high-throughput DDS helper, include the SIMD path.
#ifdef BQ_ENABLE_DDS
  #include "btrig.h"
  #include "btrig_simd.h"
#endif

namespace bq {

// ----------------------------------------------------------------------------
// Constants
// ----------------------------------------------------------------------------
constexpr long double TAU_L = 6.283185307179586476925286766559005768394L; // 2π

// ----------------------------------------------------------------------------
// PhaseQ: fixed-point turns in [0,1) with FRAC_BITS fractional bits
// ----------------------------------------------------------------------------
//  - Storage is unsigned; arithmetic is modulo 2^FRAC_BITS (pure permutation).
//  - Conversions use long double on host.
// ----------------------------------------------------------------------------

template<int FRAC_BITS = 64>
struct PhaseQ {
    static_assert(FRAC_BITS >= 1 && FRAC_BITS <= 128,
                  "FRAC_BITS must be in [1,128]");

#if defined(__SIZEOF_INT128__)
    using u128 = unsigned __int128;
#else
    // If you need MSVC or no-__int128 later, provide a portable u128.
    // For now we assume native support where FRAC_BITS>64 is used.
#endif
    using u64  = std::uint64_t;
    using stor = std::conditional_t<(FRAC_BITS <= 64), u64, u128>;

    stor frac = 0; // [0, 2^FRAC_BITS)

    // ---- Constructors / factories ------------------------------------------
    static constexpr PhaseQ from_turns_ld(long double turns) {
        long double ip; long double f = ::modfl(turns, &ip);
        if (f < 0.0L) f += 1.0L;
        long double scaled = ::ldexpl(f, FRAC_BITS);
        PhaseQ out{};
        if constexpr (FRAC_BITS <= 64) {
            if (scaled <= 0.0L) out.frac = 0;
            else if (scaled >= ::ldexpl(1.0L, FRAC_BITS)) out.frac = 0; // wrap
            else out.frac = (u64)(scaled + 0.5L);
        } else {
#if defined(__SIZEOF_INT128__)
            out.frac = (u128)(scaled + 0.5L);
#else
            // Not supported without native 128 here.
#endif
        }
        return out;
    }
    static constexpr PhaseQ from_radians_ld(long double rad) {
        return from_turns_ld(rad / TAU_L);
    }
    static constexpr PhaseQ from_degrees_ld(long double deg) {
        return from_turns_ld(deg / 360.0L);
    }

    // ---- Queries ------------------------------------------------------------
    constexpr PhaseQ half() const {
        // +1/2 turn (0.5) in Q(FRAC_BITS)
        if constexpr (FRAC_BITS <= 64) {
            return add_frac((stor) (u64(1) << (FRAC_BITS - 1)));
        } else {
#if defined(__SIZEOF_INT128__)
            return add_frac( (stor) ((u128)1 << (FRAC_BITS - 1)) );
#else
            return *this; // noop if 128 unavailable
#endif
        }
    }

    long double to_turns_ld() const {
        if constexpr (FRAC_BITS <= 64) {
            return (long double)((long double)( (u64)frac ) * std::ldexpl(1.0L, -FRAC_BITS));
        } else {
#if defined(__SIZEOF_INT128__)
            u64 hi = (u64)(((u128)frac) >> 64);
            u64 lo = (u64)((u64)frac);
            long double val = (long double)hi * std::ldexpl(1.0L, 64) + (long double)lo;
            return val * std::ldexpl(1.0L, -FRAC_BITS);
#else
            return 0.0L;
#endif
        }
    }
    long double to_radians_ld() const { return to_turns_ld() * TAU_L; }

    // ---- Modular ops --------------------------------------------------------
    constexpr PhaseQ add(const PhaseQ& d) const { PhaseQ o{}; o.frac = frac + d.frac; return o; }
    constexpr PhaseQ sub(const PhaseQ& d) const { PhaseQ o{}; o.frac = frac - d.frac; return o; }
    constexpr PhaseQ add_frac(stor k) const { PhaseQ o{}; o.frac = frac + k; return o; }
};

// Friendly aliases
using Phase32 = PhaseQ<32>;
using Phase64 = PhaseQ<64>;

// Quick helpers
inline Phase64 phase_from_turns(long double t) { return Phase64::from_turns_ld(t); }
inline Phase64 phase_from_degrees(long double d) { return Phase64::from_degrees_ld(d); }
inline Phase64 phase_from_radians(long double r) { return Phase64::from_radians_ld(r); }

// ----------------------------------------------------------------------------
// OpenQASM 3 Emitter (minimal, but now with raw-expression forms + register name)
// ----------------------------------------------------------------------------
class OpenQASM3Emitter {
public:
    explicit OpenQASM3Emitter(std::string default_reg_name = "q")
    : reg_name_(std::move(default_reg_name)) { header(); }

    void header() {
        lines_.clear();
        lines_.push_back("OPENQASM 3.0;");
        lines_.push_back("include \"stdgates.inc\";");
    }

    // Declare a qubit array; sets the active register name for subsequent ops.
    std::vector<int> qreg(const std::string& name, int n) {
        reg_name_ = name;
        lines_.push_back("qubit[" + std::to_string(n) + "] " + name + ";");
        last_size_ = n;
        std::vector<int> out(n); for (int i=0;i<n;++i) out[i]=i; return out;
    }

    // Single-qubit rotations (numeric radian)
    void rz(int q, long double a) { line("rz(" + f(a) + ") " + qref(q) + "; "); }
    void rx(int q, long double a) { line("rx(" + f(a) + ") " + qref(q) + "; "); }
    void ry(int q, long double a) { line("ry(" + f(a) + ") " + qref(q) + "; "); }

    // Raw-expression rotations (exact dyadics etc.)
    void rz_expr(int q, const std::string& expr) { line("rz(" + expr + ") " + qref(q) + "; "); }
    void rx_expr(int q, const std::string& expr) { line("rx(" + expr + ") " + qref(q) + "; "); }
    void ry_expr(int q, const std::string& expr) { line("ry(" + expr + ") " + qref(q) + "; "); }

    // Common gates
    void x (int q)              { line("x "  + qref(q) + "; "); }
    void cx(int c, int t)       { line("cx " + qref(c) + ", " + qref(t) + "; "); }
    void ccx(int c0,int c1,int t){ line("ctrl(2) @ x " + qref(c0) + ", " + qref(c1) + ", " + qref(t) + "; "); }
    void crz(int c, int t, long double a){ line("ctrl @ rz(" + f(a) + ") " + qref(c) + ", " + qref(t) + "; "); }

    void barrier() { line("barrier;"); }
    void comment(const std::string& s){ line("// " + s); }

    const char* c_str() { buffer_ = str(); return buffer_.c_str(); }
    std::string str() const { std::string out; out.reserve(1024); for (auto& s:lines_){ out+=s; out+='\n'; } return out; }
    void clear() { lines_.clear(); }

private:
    std::vector<std::string> lines_;
    std::string reg_name_ = "q";
    int last_size_ = 0;
    mutable std::string buffer_;

    std::string qref(int idx) const { return reg_name_ + "[" + std::to_string(idx) + "]"; }
    static std::string f(long double x){ std::ostringstream oss; oss.setf(std::ios::scientific); oss.precision(18); oss << (double)x; return oss.str(); }
    void line(const std::string& s){ lines_.push_back(s); }
};

// ----------------------------------------------------------------------------
// Angle emission modes
// ----------------------------------------------------------------------------
enum class AngleEmit { Float64, Dyadic };

namespace detail {
    // Denominator 2^FB as string (supports FB<=64 exactly)
    template<int FB>
    inline std::string pow2_str() {
        static_assert(FB>=1 && FB<=64, "Dyadic emitter supports FB<=64");
        if constexpr (FB < 64) return std::to_string( (std::uint64_t)1 << FB );
        else return std::string("18446744073709551616"); // 2^64
    }

    template<int FB>
    inline std::string rz_dyadic_expr(std::uint64_t k) {
        return std::string("(2*pi*") + std::to_string(k) + ")/" + pow2_str<FB>();
    }
}

// Exact dyadic emission (no rounding)
// Note: supports FRAC_BITS <= 64 for now.

template<int FB>
inline void emit_rz_dyadic(OpenQASM3Emitter& em, int q, const PhaseQ<FB>& p) {
    static_assert(FB <= 64, "emit_rz_dyadic supports FB<=64");
    em.rz_expr(q, detail::rz_dyadic_expr<FB>((std::uint64_t)p.frac));
}

template<int FB>
inline void emit_rx_dyadic(OpenQASM3Emitter& em, int q, const PhaseQ<FB>& p) {
    static_assert(FB <= 64, "emit_rx_dyadic supports FB<=64");
    em.rx_expr(q, detail::rz_dyadic_expr<FB>((std::uint64_t)p.frac));
}

template<int FB>
inline void emit_ry_dyadic(OpenQASM3Emitter& em, int q, const PhaseQ<FB>& p) {
    static_assert(FB <= 64, "emit_ry_dyadic supports FB<=64");
    em.ry_expr(q, detail::rz_dyadic_expr<FB>((std::uint64_t)p.frac));
}

// Numeric emission (host long double → rad → float literal)

template<int FB>
inline void emit_rz(OpenQASM3Emitter& em, int q, const PhaseQ<FB>& p, AngleEmit mode = AngleEmit::Float64) {
    if (mode == AngleEmit::Dyadic) emit_rz_dyadic(em, q, p); else em.rz(q, p.to_radians_ld());
}

template<int FB>
inline void emit_rx(OpenQASM3Emitter& em, int q, const PhaseQ<FB>& p, AngleEmit mode = AngleEmit::Float64) {
    if (mode == AngleEmit::Dyadic) emit_rx_dyadic(em, q, p); else em.rx(q, p.to_radians_ld());
}

template<int FB>
inline void emit_ry(OpenQASM3Emitter& em, int q, const PhaseQ<FB>& p, AngleEmit mode = AngleEmit::Float64) {
    if (mode == AngleEmit::Dyadic) emit_ry_dyadic(em, q, p); else em.ry(q, p.to_radians_ld());
}

// Vector helper

template<int FB>
inline void emit_rz_list(OpenQASM3Emitter& em,
                         const std::vector<int>& qidx,
                         const std::vector< PhaseQ<FB> >& th,
                         AngleEmit mode = AngleEmit::Float64) {
    size_t m = std::min(qidx.size(), th.size());
    for (size_t k=0;k<m;++k) emit_rz(em, qidx[k], th[k], mode);
}

// ----------------------------------------------------------------------------
// Phase schedules
// ----------------------------------------------------------------------------
// QFT-style: [1/2, 1/4, 1/8, ..., 1/2^nbits]

template<int FB>
inline std::vector< PhaseQ<FB> > make_phase_qft(int nbits) {
    std::vector< PhaseQ<FB> > out; out.reserve(nbits);
    for (int j=1;j<=nbits;++j) {
        long double turns = std::ldexpl(1.0L, -j); // 1/2^j
        out.push_back( PhaseQ<FB>::from_turns_ld(turns) );
    }
    return out;
}

// Linear phases in turns (already had a version; keep here for completeness)

template<int FB>
inline std::vector< PhaseQ<FB> > make_phase_linear(std::size_t count,
                                                   long double turns_start,
                                                   long double turns_step) {
    std::vector< PhaseQ<FB> > out; out.reserve(count);
    long double t = turns_start;
    for (std::size_t i = 0; i < count; ++i) { out.push_back( PhaseQ<FB>::from_turns_ld(t) ); t += turns_step; }
    return out;
}

// ----------------------------------------------------------------------------
// (Optional) DDS streamer — host-side high-throughput sincos using btrig_simd
// ----------------------------------------------------------------------------
#ifdef BQ_ENABLE_DDS
inline void dds_stream_q64(btrig::uq64 phase0,
                           btrig::uq64 dphase,
                           std::size_t N,
                           double* cosOut,
                           double* sinOut,
                           bool precise = false) {
    btrig::uq64 p = phase0;
    constexpr std::size_t CH = 2048; // chunk size for cache friendliness
    std::vector<btrig::uq64> buf; buf.resize(std::min(CH, N));
    std::size_t i = 0;
    while (i < N) {
        std::size_t m = std::min(CH, N - i);
        for (std::size_t k = 0; k < m; ++k) { buf[k] = p; p += dphase; }
        btrig::simd::sincos_q64(buf.data(), (int)m, precise, cosOut + i, sinOut + i);
        i += m;
    }
}
#endif

// ----------------------------------------------------------------------------
// Notes
// ----------------------------------------------------------------------------
// • Use AngleEmit::Dyadic for exact textual angles in OpenQASM 3 outputs—great
//   for papers/validation. Float64 is fine for most workflows and is readable.
// • Keep heavy math on the host (btrig/btrig_simd). On the QPU side, prefer
//   parameterized unitaries over trying to compute trig inside a circuit.
// ----------------------------------------------------------------------------

} // namespace bq

#endif // BTRIG_QUANTUM_H
