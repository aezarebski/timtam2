package timtam;

import beast.base.util.Binomial;

/**
 * <p>Representation of the number of hidden lineages used by {@link TimTam}.
 * Originally this was strictly a negative binomial distribution but to
 * accommodate additional features has been extended to include degenerate
 * distributions where the probability mass is concentrated on a single
 * value.</p>
 *
 * @author Alexander E. Zarebski
 */
public class HiddenLineageDist {

    // We track the logarithm of the mean and variance because these can vary
    // substantially.
    private double lnMean;
    private double lnVariance;

    // When using the negative binomial the r and p parameterization is
    // convenient, so we keep the distribution in terms of these too.
    private double lnP;
    private double ln1mP;
    private double lnR;

    // The following are useful as a way to check what kind of distribution the
    // object is representing.
    private boolean isZero;
    private boolean isDegenerate;

    /**
     * It is not entirely clear why we need to have an empty constructor, but it
     * seems to be required for something relating to introspection.
     */
    public HiddenLineageDist() {
    }

    /**
     * The logarithm of the probability generating function using the supplied
     * parameters.
     *
     * @param z variable of the generating function
     * @param r parameter of the negative distribution
     * @param p parameter of the distribution
     *
     * @return the logarithm of the probability generating function.
     */
    public double lnPGF(double z, double r, double p) {
        double ln1mp = Math.log(1 - p);
        double ln1mpz = Math.log(1 - p * z);
        return r * (ln1mp - ln1mpz);
    }

    /**
     * The logarithm of the probability generating function for this instance of
     * the negative binomial distribution.
     *
     * @param z variable of the generating function
     *
     * @return the logarithm of the probability generating function.
     *
     * @throws RuntimeException if the distribution is degenerate.
     */
    public double lnPGF(double z) {
        if (!this.isDegenerate) {
            double p = Math.exp(this.lnP);
            double r = Math.exp(this.lnR);
            return lnPGF(z, r, p);
        } else if (this.isZero) {
            return 0.0;
        }   else {
            return Math.exp(this.lnMean) * Math.log(z);
        }
    }

    /**
     * The logarithm of the n-th derivative of the probability generating
     * function using the supplied parameters.
     *
     * @param n the number of partial derivatives
     * @param z variable of the generating function
     * @param r parameter of the negative distribution
     * @param p parameter of the distribution
     *
     * @return the logarithm of the n-th partial derivative of the probability
     * generating function
     */
    public double lnPGFDash(int n, double z, double r, double p) {
        if (n > 0) {
            double lnP = Math.log(p);
            double ln1mP = Math.log(1 - p);
            return Numerics.lnPochhammer(r, n) + n * (lnP - ln1mP) + lnPGF(z, r + n, p);
        } else if (n == 0) {
            return lnPGF(z, r, p);
        } else {
            throw new IllegalArgumentException("lnPGFDash expected n >= 0 but got " + n);
        }
    }

    /**
     * The logarithm of the n-th derivative of the probability generating
     * function for this instance of the negative binomial distribution.
     *
     * @param n the number of partial derivatives
     * @param z variable of the generating function
     *
     * @return the logarithm of the n-th partial derivative of the probability
     * generating function
     */
    public double lnPGFDash(int n, double z) throws IllegalArgumentException {
        if (!this.isDegenerate) {
            double p = Math.exp(this.lnP);
            double r = Math.exp(this.lnR);
            return lnPGFDash(n, z, r, p);
        } else if (this.isZero) {
            if (n > 0) {
                return Double.NEGATIVE_INFINITY;
            } else if (n == 0) {
                return 0.0;
            } else {
                throw new IllegalArgumentException("lnPGFDash expected n >= 0 but got " + n);
            }
        } else {
            if (n > 0) {
                double m = Math.exp(this.lnMean);
                return Numerics.lnPochhammer(m,n) + (m-n)*Math.log(z);
            } else if (n == 0) {
                return Math.exp(this.lnMean) * Math.log(z);
            } else {
                throw new IllegalArgumentException("lnPGFDash expected n >= 0 but got " + n);
            }
        }
    }

    public double lnPGFDash1(double z) {
        return lnPGFDash(1, z);
    }

    public double lnPGFDash2(double z) {
        return lnPGFDash(2, z);
    }


    public double getLnMean() {
        return lnMean;
    }

    public double getLnVariance() {
        return lnVariance;
    }


    /**
     * Update the state of the distribution of hidden lineages based on the mean
     * and variance.
     *
     * <b>NOTE:</b> The current implementation does not support the case where the
     * variance is positive and less than the mean because this situation cannot
     * be represented by either a negative binomial distribution or a degenerate
     * distribution.
     *
     * @param lnMean the natural logarithm of the mean of the number of hidden
     * lineages.
     * @param lnVariance the natural logarithm of the variance in the number of
     * hidden lineages.
     * @throws RuntimeException
     */
    public void setLnMeanAndLnVariance(double lnMean, double lnVariance) throws RuntimeException {
        this.lnMean = lnMean;
        this.lnVariance = lnVariance;

        if (Double.isFinite(this.lnVariance) & this.lnMean >= this.lnVariance) {
            throw new RuntimeException(
                "HiddenLineageDist does not support non-degenerate distributions where the " +
                "variance is less than the mean. You may be able to fix this by adjusting the " +
                "times at which you estimate the prevalence."
            );
        }

        double lnvmm = Math.log(Math.exp(lnVariance) - Math.exp(lnMean));
        this.lnP = lnvmm - lnVariance;
        this.ln1mP = Math.log(1 - Math.exp(lnP));
        this.lnR = 2 * lnMean - lnvmm;

        this.isZero = false;
        this.isDegenerate = false;
    }

    public void setLnPAndLnR(double lnP, double lnR) {
        this.lnP = lnP;
        this.ln1mP = Math.log(1 - Math.exp(lnP));
        this.lnR = lnR;

        // TODO Work out a safer way to compute this.
        this.lnMean = lnP + lnR - this.ln1mP;
        this.lnVariance = this.lnMean - this.ln1mP;

        this.isZero = false;
        this.isDegenerate = false;
    }

    public boolean getIsZero() {
        return this.isZero;
    }

    public double getLnP() {
        return this.lnP;
    }

    public double getLnR() {
        return this.lnR;
    }

    /**
     * Return the log-probability of this distribution.
     *
     * @param n the argument to the PMF
     * @return the log-probability of this distribution.
     */
    public double lnPMF(int n) {
        if (!this.isDegenerate) {
            // The distribution must follow a negative binomial distribution.
            double r = Math.exp(this.lnR);
            return Numerics.lnChoose(n + r - 1, n) +
                r * this.ln1mP +
                n * this.lnP;
        } else if (this.isZero) {
            // The distribution must be degenerate at zero.
            return (n == 0) ?
                0.0 :
                Double.NEGATIVE_INFINITY;
        } else {
            return (n == Math.round(Math.exp(this.lnMean))) ?
                0.0 :
                Double.NEGATIVE_INFINITY;
        }
    }

    /**
     * Set this object to represent a degenerate distribution with a point mass
     * at the given integer value.
     *
     * @param x the location of the point mass
     */
    public void setIsDegenerate(int x) {
        if (x > 0) {
            this.lnMean = Math.log(x);
            this.isZero = false;
        } else if (x == 0) {
            this.lnMean = Double.NEGATIVE_INFINITY;
            this.isZero = true;
        } else {
            throw new RuntimeException("setIsDegenerate expects a non-negative integer but got " + x);
        }
        this.lnVariance = Double.NEGATIVE_INFINITY;
        this.lnP = Double.NaN;
        this.ln1mP = Double.NaN;
        this.lnR = Double.NaN;
        this.isDegenerate = true;
    }

    public boolean getIsDegenerate() {
        return this.isDegenerate;
    }
}
