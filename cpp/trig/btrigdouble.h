namespace btrigdouble {

    // ---- Packed per-anchor with split θ ----
    struct Anchor { double thHi, thLo, c, s; };
    alignas(64) static Anchor gAnchors[btrig::NUM_ANCHORS_FULL];

    // Split TAU (used for q*TAU reduction without 1/TAU multiplies)
    static double TAU_HI_S = 0.0, TAU_LO_S = 0.0, TAU_SUM_S = 0.0;

    // Explicit one-time init (call this once after btrig::initFullCircleAnchors())
    static inline void init() {
        static int inited = 0; if (inited) return;
        const double C = 134217729.0; // 2^27+1

        // split(TAU)
        double t_tau = C * btrig::TAU;
        TAU_HI_S = t_tau - (t_tau - btrig::TAU);
        TAU_LO_S = btrig::TAU - TAU_HI_S;
        TAU_SUM_S = TAU_HI_S + TAU_LO_S;

        // pack anchors with split θ
        for (int j = 0; j < btrig::NUM_ANCHORS_FULL; ++j) {
            double th = btrig::fullCircleAnchors[j][0];
            double t = C * th;
            gAnchors[j].thHi = t - (t - th);
            gAnchors[j].thLo = th - gAnchors[j].thHi;
            gAnchors[j].c = btrig::fullCircleAnchors[j][1];
            gAnchors[j].s = btrig::fullCircleAnchors[j][2];
        }
        inited = 1;
    }

    // Precise-only sincos with robust indexing
    static __forceinline void sincos(double angle, bool /*precise*/,
        double& cosOut, double& sinOut)
    {
        // --- robust index: double-double tf = angle * ONE_OVER_STEP ---
        double tf_hi = angle * btrig::ONE_OVER_STEP;
        double tf_lo = fma(angle, btrig::ONE_OVER_STEP, -tf_hi); // exact residual

        // fast floor of tf_hi
        int i = (int)tf_hi;
        i -= ((tf_hi < 0.0) & ((double)i != tf_hi));

        // use lo-part to fix off-by-1 at seams
        double frac = (tf_hi - (double)i) + tf_lo;   // (-eps, 1+eps)
        if (frac < 0.0) { --i; /*frac += 1.0;*/ }
        else if (frac >= 1.0) { ++i; /*frac -= 1.0;*/ }

        // turn count and reduced angle via split TAU (no 1/TAU multiplies)
        constexpr int N = btrig::NUM_ANCHORS_FULL;
        constexpr int MSK = btrig::ANCHOR_MASK;

        const int i_mod = (i & MSK);           // 0..N-1 (exact for power-of-two N)
        const int q = (i - i_mod) / N;     // exact integer division

        double ang = fma(-(double)q, TAU_HI_S, angle);
        ang = fma(-(double)q, TAU_LO_S, ang);   // now 0 <= ang < TAU

        // seam-cancel pivot: previous cell
        const int idx = (i_mod - 1) & MSK;
        const Anchor& A = gAnchors[idx];

        // residual about split-θ pivot
        double d = fma(-1.0, A.thHi, ang);
        d = fma(-1.0, A.thLo, d);

        // rare wrap fix in first cell
        if (i_mod == 0 && d < 0.0) d += TAU_SUM_S;

        // C8/S7 micro
        const double d2 = d * d;

        const double C2 = -1.0 / 2.0, C4 = 1.0 / 24.0, C6 = -1.0 / 720.0, C8 = 1.0 / 40320.0;
        double pc = fma(d2, C8, C6);
        pc = fma(d2, pc, C4);
        pc = fma(d2, pc, C2);
        const double cosd = fma(d2, pc, 1.0);

        const double S3 = -1.0 / 6.0, S5 = 1.0 / 120.0, S7 = -1.0 / 5040.0;
        double ps = fma(d2, S7, S5);
        ps = fma(d2, ps, S3);
        const double sind = d * fma(d2, ps, 1.0);

        // rotate from anchor
        cosOut = fma(A.c, cosd, -(A.s * sind));
        sinOut = fma(A.s, cosd, (A.c * sind));
    }

    // Keep API parity and ensure precise-only path is always used
    static __forceinline void sin(double angle, bool /*precise*/, double& outSin) {
        double c; sincos(angle, true, c, outSin);
    }
    static __forceinline void cos(double angle, bool /*precise*/, double& outCos) {
        double s; sincos(angle, true, outCos, s);
    }
}
