package org.brennon.com.server;

import static org.brennon.com.server.TaylorTrig.NUM_ANCHORS_QUADRANT;
import static org.brennon.com.server.TaylorTrig.TAU;

public class TrigBenchmark {

    static boolean precise = false;

    // Main function to test the code
    public static void main(String[] args) {
        // Run the accuracy benchmark
        benchmarkSpeed(1); // No depth needed here anymore
        benchmarkAccuracy();
    }

    // Speed Benchmark (Taylor Series vs Trig)
    public static void benchmarkSpeed(int numSamples) {
        double[] angles = new double[numSamples];
        double trigSum = 0;
        double taylorSum = 0;

        // Prepare angles
        for (int i = 0; i < numSamples; i++) {
            angles[i] = Math.random() * TAU;
        }

        // Warmup for Math.cos/sin
        for (int i = 0; i < 1000000; i++) {  // Run a warmup loop
            double angle = angles[i % numSamples];
            Math.cos(angle);
            Math.sin(angle);
        }

        // Warmup for Taylor-based rotation
        for (int i = 0; i < 1000000; i++) {  // Run a warmup loop
            double vTaylor[] = TaylorTrig.getUnitVectorFromAngle(angles[i % numSamples], precise);  // Using precise mode
        }

        // Timing Math.cos/sin (sum of components)
        long start = System.nanoTime();
        for (int i = 0; i < numSamples; i++) {
            double angle = angles[i];
            double xTrig = Math.cos(angle);
            double yTrig = Math.sin(angle);
            trigSum += xTrig + yTrig;  // Sum of components to prevent optimization
        }
        long end = System.nanoTime();
        double trigTime = (end - start) / 1e6; // Convert to milliseconds

        // Timing Taylor-based rotation (sum of components)
        start = System.nanoTime();
        for (int i = 0; i < numSamples; i++) {
            double vTaylor[] = TaylorTrig.getUnitVectorFromAngle(angles[i], precise);  // Using precise mode
            taylorSum += vTaylor[0] + vTaylor[1]; // Sum of components
        }
        end = System.nanoTime();
        double taylorTime = (end - start) / 1e6; // Convert to milliseconds

        // Log the results with individual speeds
        System.out.printf("--- SPEED TEST (%d samples) ---\n", numSamples);  // No depth in the output anymore
        System.out.printf("Native Math.cos/sin (sum of x + y):  %.2f ms\n", trigTime);
        System.out.printf("Custom Taylor rotation (sum of x + y): %.2f ms\n", taylorTime);
        System.out.printf("Trig checksum (sum of x + y):        %.12f\n", trigSum);
        System.out.printf("Taylor checksum (sum of x + y):      %.12f\n", taylorSum);
    }

    // Benchmark Accuracy (Taylor vs Math.cos/sin)
    public static void benchmarkAccuracy() {
        double[] testAngles = {0.575, 1.111, 2.789, 3.333, 4.444, 5.555};

        System.out.println("--- ACCURACY TEST vs Math.cos/sin (NUM_ANCHORS_QUADRANT=" + NUM_ANCHORS_QUADRANT + ") ---");

        // Warmup for Math.cos/sin
        for (double angle : testAngles) {
            Math.cos(angle);
            Math.sin(angle);
        }

        // Warmup for Taylor-based rotation
        for (double angle : testAngles) {
            double[] approx1 = TaylorTrig.getUnitVectorFromAngle(angle, false);  // Inline Taylor expansion
        }

        for (double angle : testAngles) {
            // Normalize the angle to [0, 2π]
            double normalizedAngle = (angle % TAU + TAU) % TAU;

            // Get approximations from Taylor method
            double[] approx1 = TaylorTrig.getUnitVectorFromAngle(normalizedAngle, false);  // Inline Taylor expansion

            // Get exact values from Math.cos and Math.sin
            double exactX = Math.cos(normalizedAngle);
            double exactY = Math.sin(normalizedAngle);

            // Calculate errors for Taylor method
            double err1X = Math.abs(approx1[0] - exactX);
            double err1Y = Math.abs(approx1[1] - exactY);
            double err1 = Math.hypot(err1X, err1Y);  // Euclidean distance (error)

            // Log the results
            System.out.printf("angle=%.3f → Taylor err=%.2e%n", normalizedAngle, err1);
        }
    }
}
