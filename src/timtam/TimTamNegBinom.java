package timtam;

public class TimTamNegBinom {

    private double lnMean;
    private double lnVariance;
    private double lnP;
    private double ln1mP;
    private double lnR;
    private boolean isZero;

    public TimTamNegBinom() {
    }

    /**
     * The logarithm of the probability generating function using the supplied parameters.
     *
     * @param z variable of the generating function
     * @param r parameter of the negative distribution
     * @param p parameter of the distribution
     * @return the logarithm of the probability generating function.
     */
    public double lnPGF(double z, double r, double p) {
        double ln1mp = Math.log(1 - p);
        double ln1mpz = Math.log(1 - p * z);
        return r * (ln1mp - ln1mpz);
    }

    /**
     * The logarithm of the probability generating function for this instance of the negative binomial distribution.
     *
     * @param z variable of the generating function
     * @return the logarithm of the probability generating function.
     */
    public double lnPGF(double z) {
        double p = Math.exp(this.lnP);
        double r = Math.exp(this.lnR);
        return lnPGF(z, r, p);
    }

    /**
     * The logarithm of the n-th derivative of the probability generating function using the supplied parameters.
     *
     * @param n the number of partial derivatives
     * @param z variable of the generating function
     * @param r parameter of the negative distribution
     * @param p parameter of the distribution
     * @return the logarithm of the n-th partial derivative of the probability generating function
     */
    public double lnPGFDash(int n, double z, double r, double p) {
        if (n == 0) {
            return lnPGF(z, r, p);
        }
        if (n < 0) {
            throw new IllegalArgumentException("lnPGFDash expected n > 0 but got " + n);
        }

        double lnP = Math.log(p);
        double ln1mP = Math.log(1 - p);
        return lnPochhammer(r, n) + n * (lnP - ln1mP) + lnPGF(z, r + n, p);
    }

    /**
     * The logarithm of the n-th derivative of the probability generating function for this instance of the negative binomial distribution.
     *
     * @param n the number of partial derivatives
     * @param z variable of the generating function
     * @return the logarithm of the n-th partial derivative of the probability generating function
     */
    public double lnPGFDash(int n, double z) {
        double p = Math.exp(this.lnP);
        double r = Math.exp(this.lnR);
        return lnPGFDash(n, z, r, p);
    }

    public double lnPGFDash1(double z) {
        return lnPGFDash(1, z);
    }

    public double lnPGFDash2(double z) {
        return lnPGFDash(2, z);
    }

    /**
     * @return the logarithm of the Pochhammer function
     */
    public double lnPochhammer(double a, int i) {
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

    public double getLnMean() {
        return lnMean;
    }

    public double getLnVariance() {
        return lnVariance;
    }

    public void setLnMeanAndLnVariance(double lnMean, double lnVariance) {
        this.lnMean = lnMean;
        this.lnVariance = lnVariance;

        double lnvmm = Math.log(Math.exp(lnVariance) - Math.exp(lnMean));
        this.lnP = lnvmm - lnVariance;
        this.ln1mP = Math.log(1 - Math.exp(lnP));
        this.lnR = 2 * lnMean - lnvmm;

        this.isZero = false;
    }

    public void setLnPAndLnR(double lnP, double lnR) {
        this.lnP = lnP;
        this.ln1mP = Math.log(1 - Math.exp(lnP));
        this.lnR = lnR;

        // TODO Work out a safer way to compute this.
        this.lnMean = lnP + lnR - this.ln1mP;
        this.lnVariance = this.lnMean - this.ln1mP;

        this.isZero = false;
    }

    public boolean getIsZero() {
        return this.isZero;
    }

    public void setIsZero(boolean isZero) {
        if (isZero) {
            this.lnMean = Double.NEGATIVE_INFINITY;
            this.lnVariance = Double.NEGATIVE_INFINITY;
            this.lnP = Double.NaN;
            this.ln1mP = Double.NaN;
            this.lnR = Double.NaN;
            this.isZero = isZero;
        } else {
            throw new RuntimeException("setIsZero should not be used to set this.");
        }
    }

    public double getLnP() {
        return lnP;
    }

    public double getLnR() {
        return lnR;
    }

}