package timtam;

public class TimTamNegBinom {

    double lnMean;
    double lnVariance;
    double lnP;
    double ln1mP;
    double lnR;

    public void setZero() {
        isZero = true;
    }

    boolean isZero;

    public TimTamNegBinom() {
    }

    public String toString() {
        return "NegBinomial(" + Math.exp(this.lnR) + ", " + Math.exp(lnP) + ")";
    }

    /**
     *
     * @param z
     * @param r
     * @param p
     * @return the logarithm of the probability generating function.
     */
    public double lnPGF(double z, double r, double p) {
        double ln1mp = Math.log(1 - p);
        double ln1mpz = Math.log(1 - p * z);
        return r * (ln1mp - ln1mpz);
    }

    public double lnPGF(double z) {
        double p = Math.exp(this.lnP);
        double r = Math.exp(this.lnR);
        return lnPGF(z, r, p);
    }

    /**
     * @param n the number of partial derivatives
     * @param z
     * @return the logarithm of the n-th partial derivative of the probability generating function
     */
    public double lnPGFDash(int n, double z, double r, double p) {
        if (n == 0) {
            return lnPGF(z, r, p);
        }
        if (n < 0) {
            throw new IllegalArgumentException("lnPGFDash expected n > 0 but got " + n);
        }

        if (!this.isZero) {
            double lnP = Math.log(p);
            double ln1mP = Math.log(1 - p);
            return lnPochhammer(r, n) + n * (lnP - ln1mP) + lnPGF(z, r + n, p);
        } else {
            return Double.NEGATIVE_INFINITY;
        }
    }

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
     * @param a
     * @param i
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

    public TimTamNegBinom(double mean, double variance) {
        if (mean > 0 && variance > mean) {
            this.isZero = false;
            setLnMeanAndLnVariance(Math.log(mean), Math.log(variance));
        }

        if (mean == 0 && variance == 0) {
            this.isZero = true;
            this.lnMean = Double.NEGATIVE_INFINITY;
            this.lnVariance = Double.NEGATIVE_INFINITY;
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
        updateLnPAndLnR();
        if (lnMean > 0) this.isZero = false;
    }

    public double getLnP() {
        return lnP;
    }

    public double getLnR() {
        return lnR;
    }

    public void setLnPAndLnR(double lnP, double lnR) {
        this.lnP = lnP;
        this.ln1mP = Math.log(1 - Math.exp(lnP));
        this.lnR = lnR;
        updateLnMeanAndLnVariance();
        if (this.lnMean > 0) this.isZero = false;
    }

    /**
     * Update the values of p and r based on the mean and variance.
     */
    private void updateLnPAndLnR() {
        // TODO Work out a safer way to compute this.
        double lnvmm = Math.log(Math.exp(this.lnVariance) - Math.exp(this.lnMean));
        this.lnP = lnvmm - this.lnVariance;
        this.ln1mP = Math.log(1 - Math.exp(lnP));
        this.lnR = 2 * this.lnMean - lnvmm;
    }

    /**
     * Update the values of the mean and variance based on p and r.
     */
    private void updateLnMeanAndLnVariance() {
        // TODO Work out a safer way to compute this.
        this.lnMean = this.lnP + this.lnR - this.ln1mP;
        this.lnVariance = this.lnMean - this.ln1mP;
    }

}