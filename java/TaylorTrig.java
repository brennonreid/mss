// Constants
package org.brennon.com.server;

// Quadrant Transforms (functional interface)
class TaylorTrig {

    public static final double TAU = 2 * Math.PI;
    public static final double QUAD_TAU = TAU / 4;
    public static final int NUM_ANCHORS_QUADRANT = 64;
    public static final double DEG_STEP_QUADRANT = QUAD_TAU / NUM_ANCHORS_QUADRANT;
    public static final int NUM_ANCHORS_FULL = 4 * NUM_ANCHORS_QUADRANT;
    public static final double DEG_STEP_FULL = TAU / NUM_ANCHORS_FULL;
    public static final double STEP = TAU / NUM_ANCHORS_FULL;


    public TaylorTrig(){}

    @FunctionalInterface
    public interface QuadrantTransform {
        double[] transform(double x, double y);
    }

    public static final QuadrantTransform[] quadrantTransforms = new QuadrantTransform[]{
            (x, y) -> new double[]{x, y},
            (x, y) -> new double[]{-y, x},
            (x, y) -> new double[]{-x, -y},
            (x, y) -> new double[]{y, -x}
    };

    // Apply quadrant transformation
    public static double[] applyQuadrantTransform(double x, double y, int quadrant) {
        return quadrantTransforms[quadrant].transform(x, y);
    }

    // Taylor Series Expansion for Cosine
    public static double taylorCosFloat(double x, int depth) {
        x = x % TAU;  // Normalize angle
        double result = 1;  // First term (1)
        double term = 1;    // Initial term
        double x2 = x * x;  // x^2 for reuse

        for (int i = 1; i <= depth; i++) {
            term *= -x2 / ((2 * i - 1) * (2 * i));  // Calculate the next term using Horner's method
            result += term;  // Add it to the result
        }
        return result;
    }

    // Taylor Series Expansion for Sine
    public static double taylorSinFloat(double x, int depth) {
        x = x % TAU;  // Normalize angle
        double result = x;  // First term (x)
        double term = x;    // Initial term
        double x2 = x * x;  // x^2 for reuse

        for (int i = 1; i <= depth; i++) {
            term *= -x2 / ((2 * i) * (2 * i + 1));  // Calculate the next term using Horner's method
            result += term;  // Add it to the result
        }
        return result;
    }
    // Canonical Anchors (for base angle steps)
    public static class Anchor {
        double theta, cos, sin;

        public Anchor(double theta, double cos, double sin) {
            this.theta = theta;
            this.cos = cos;
            this.sin = sin;
        }
    }

    public static Anchor[] canonicalAnchors = new Anchor[NUM_ANCHORS_QUADRANT];

    static {
        for (int i = 0; i < NUM_ANCHORS_QUADRANT; i++) {
            double angle = i * DEG_STEP_QUADRANT;
            canonicalAnchors[i] = new Anchor(angle, taylorCosFloat(angle, 55), taylorSinFloat(angle, 55));
        }
    }

    // Full Circle Anchors (based on canonical anchors)
    public static Anchor[] fullCircleAnchors = new Anchor[NUM_ANCHORS_FULL];

    static {
        for (int i = 0; i < NUM_ANCHORS_FULL; i++) {
            double angle = i * DEG_STEP_FULL;
            int quadrant = (int) Math.floor(angle / QUAD_TAU);
            double theta = angle - quadrant * QUAD_TAU;
            int canonicalIndex = (int) Math.round(theta / DEG_STEP_QUADRANT);
            Anchor canonicalAnchor = canonicalAnchors[Math.min(canonicalIndex, NUM_ANCHORS_QUADRANT - 1)];
            double[] transformed = applyQuadrantTransform(canonicalAnchor.cos, canonicalAnchor.sin, quadrant);
            fullCircleAnchors[i] = new Anchor(angle, transformed[0], transformed[1]);
        }
    }

    // Corrected Symbolic Inflate function with proper higher-order terms (3 and 4)
    public static double[] getUnitVectorFromAngle(double angleRadians, boolean precise) {
        double TAU = Math.PI * 2;
        // Assuming DEG_STEP_FULL and NUM_ANCHORS_FULL are defined elsewhere.
        double DEG_STEP_FULL = TAU / NUM_ANCHORS_FULL;
        double relativeAngle = angleRadians;

        relativeAngle %= TAU;
        if (relativeAngle < 0) relativeAngle += TAU;

        int zoneIndex = (int) Math.floor(relativeAngle / DEG_STEP_FULL) % NUM_ANCHORS_FULL;
        Anchor anchor = fullCircleAnchors[zoneIndex];

        double delta = (relativeAngle - anchor.theta);
        double deltaSquared = delta * delta;

        // horner term 3 for cos
        double cosRel = 1 + deltaSquared * (-0.5 + deltaSquared * (1.0 / 24.0) + deltaSquared * deltaSquared * (-1.0 / 720.0));
        if (precise) {
            // horner term 4 for cos
            cosRel += deltaSquared * deltaSquared * deltaSquared * deltaSquared * (1.0 / 40320.0);
        }

        // horner term 3 for sin
        double sinRel = delta + delta * deltaSquared * (-1.0 / 6.0 + deltaSquared * (1.0 / 120.0) + deltaSquared * deltaSquared * (-1.0 / 5040.0));
        if (precise) {
            // horner term 4 for sin
            sinRel += delta * deltaSquared * deltaSquared * deltaSquared * deltaSquared * (1.0 / 362880.0);
        }

        // Rotation formula
        double x = anchor.cos * cosRel - anchor.sin * sinRel;
        double y = anchor.cos * sinRel + anchor.sin * cosRel;
        return new double[]{x, y};
    }



}
