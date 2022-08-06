package timtam;

import beast.math.GammaFunction;

/**
 * Useful numerical functions.
 *
 * @author Alexander E. Zarebski
 */
public class Numerics {

    /**
     * The logarithm of the <a href="https://en.wikipedia.org/wiki/Falling_and_rising_factorials">Pochhammer function</a>.
     */
    public static double lnPochhammer(double a, int i) {
        if (i > 0) {
            double tmp = 0;
            for (int j = i; j > 0; j--) {
                tmp += Math.log(a + j - 1);
            }
            return tmp;
        } else if (i == 0) {
            return 0.0;
        } else {
            throw new IllegalArgumentException("logPochhammer expects i >= 0 but got " + i);
        }
    }

    /**
     * The logarithm of the binomial coefficient where the first parameter can take real values.
     */
    public static double lnChoose(final double n, final int k) {
        return GammaFunction.lnGamma(n + 1.0) - GammaFunction.lnGamma(k + 1.0)
            - GammaFunction.lnGamma(n - k + 1.0);
    }


    /**
     * The sum of the values in the array
     */
    public static double arraySum(double[] xs) {
        double tmp = 0.0;
        for (double x : xs) {
            tmp += x;
        }
        return(tmp);
    }

    /**
     * The sum of the values in the array
     */
    public static int arraySum(int[] xs) {
        int tmp = 0;
        for (int x : xs) {
            tmp += x;
        }
        return(tmp);
    }

    /**
     * LogSumExp
     *
     * <p>This function implements a numerically safe way to take the logarithm
     * of the sum of exponentials.</p>
     *
     * @see <a href="https://en.wikipedia.org/wiki/LogSumExp">Wikipedia page.</a>
     */
    public static double logSumExp(double[] xs) {
        return logSumExp(xs, xs.length);
    }

    /**
     * This is a nitty-gritty log-sum-exp which is useful when fine-tuning
     * memory usage by avoiding streams.
     *
     * @param xs array
     * @param n number of leading entries to use
     */
    public static double logSumExp(double[] xs, int n) {
        double xMax = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < n; i++) {
            if (xs[i] > xMax) {
                xMax = xs[i];
            }
        }
        double tmp = 0;
        for (int i = 0; i < n; i++) {
            tmp += Math.exp(xs[i] - xMax);
        }
        return xMax + Math.log(tmp);
    }
}
