package beast.evolution.speciation;


import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.birthdeath.EventType;
import beast.evolution.tree.birthdeath.TreeWithPointProcess;

import java.util.Arrays;

/**
 * Tree prior for birth-death-sampling while tracking the distribution of hidden
 * lineages. This used BirthDeathSerialSampling.java as a starting point. There
 * is a nested class, NegativeBinomial, which is the approximation used for the
 * number of unobserved lineages in the likelihood.
 *
 * @author Alexander Zarebski
 */
public class TimTam extends TreeDistribution {

    // birth rate
    RealParameter lambda;

    // death rate
    RealParameter mu;

    // serial sampling rate
    RealParameter psi;

    // extant sampling proportion
    RealParameter p;

    // occurrence rate
    RealParameter omega;

    // length of the edge from origin to MRCA node.
    RealParameter rootLength;

    // the times at which a scheduled sequenced sample was attempted
    TraitSet catastropheTimes;

    // the times at which there was an occurrence sample.
    TraitSet points;

    // we use this attribute to accumulate the log-likelihood in the calculation.
    private double lnL;

    // this is used to track the number of lineages in the reconstructed tree as we process the data.
    private int k;

    // this is used to track the distribution of hidden lineages as we process the data.
    private NegativeBinomial nb;

    // Unclear why it is necessary, but BEAST expects there to be a zero-argument
    // constructor and if there isn't one it freaks out.
    public TimTam() {

    }

    public TimTam(
            RealParameter lambda,
            RealParameter mu,
            RealParameter psi,
            RealParameter p,
            RealParameter omega,
            RealParameter rootLength,
            TraitSet catastropheTimes,
            TraitSet points) {
        this("timTamModel", lambda, mu, psi, p, omega, rootLength, catastropheTimes, points);
    }

    final public Input<RealParameter> lambdaInput = new Input<>("lambda", "the birth rate of new infections");
    final public Input<RealParameter> muInput = new Input<>("mu", "the death rate");
    final public Input<RealParameter> psiInput = new Input<>("psi", "the sampling rate");
    final public Input<RealParameter> pInput = new Input<>("p", "the probability of sampling extant lineages");
    final public Input<RealParameter> omegaInput = new Input<>("omega", "the occurrence rate");
    final public Input<RealParameter> rootLengthInput = new Input<>("rootLength", "the length of the edge between the origin and the MRCA");
    final public Input<TraitSet> catastropheTimesInput = new Input<>("catastropheTimes", "the times at which a scheduled sequenced sample was attempted");
    final public Input<TraitSet> pointsInput = new Input<>("points", "the points in the point process");

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        this.lambda = lambdaInput.get();
        lambda.setBounds(0.0, Double.POSITIVE_INFINITY);

        this.mu = muInput.get();
        mu.setBounds(0.0, Double.POSITIVE_INFINITY);

        this.psi = psiInput.get();
        psi.setBounds(0.0, Double.POSITIVE_INFINITY);

        this.p = pInput.get();
        p.setBounds(0.0, 1.0);

        this.omega = omegaInput.get();
        omega.setBounds(0.0, Double.POSITIVE_INFINITY);

        this.rootLength= rootLengthInput.get();
        rootLength.setBounds(0.0, Double.POSITIVE_INFINITY);

        this.catastropheTimes = catastropheTimesInput.get();

        this.points = pointsInput.get();
    }

    public TimTam(
            String modelName,
            RealParameter lambda,
            RealParameter mu,
            RealParameter psi,
            RealParameter p,
            RealParameter omega,
            RealParameter rootLength,
            TraitSet catastropheTimes,
            TraitSet points) {

        this.lambda = lambda;
        lambda.setBounds(0.0, Double.POSITIVE_INFINITY);

        this.mu = mu;
        mu.setBounds(0.0, Double.POSITIVE_INFINITY);

        this.psi = psi;
        psi.setBounds(0.0, Double.POSITIVE_INFINITY);

        this.p = p;
        p.setBounds(0.0, 1.0);

        this.omega = omega;
        omega.setBounds(0.0, Double.POSITIVE_INFINITY);

        this.rootLength= rootLength;
        rootLength.setBounds(0.0, Double.POSITIVE_INFINITY);

        this.nb = new NegativeBinomial();

        this.catastropheTimes = catastropheTimes;

        this.points = points;
    }

    public double birth() {
        return lambda.getValue(0);
    }

    public double death() {
        return mu.getValue(0);
    }

    public double psi() {
        return psi.getValue(0);
    }

    public double omega() {
        return omega.getValue(0);
    }

    public double p() {
        return p.getValue(0);
    }


    @Override
    public double calculateLogP() {
        calculateTreeLogLikelihood((Tree) treeInput.get());
        return this.lnL;
    }

    /**
     * Generic likelihood calculation
     *
     * @param tree the tree to calculate likelihood of
     * @return log-likelihood of density
     */
    public final void calculateTreeLogLikelihood(Tree tree) {

        TreeWithPointProcess ti = new TreeWithPointProcess(rootLength, tree, points, catastropheTimes);
        int numIntervals = ti.getIntervalCount();

        this.lnL = 0.0;
        // TODO this assumes a fixed initial condition which should be movied
        //  into something that gets specified by the client.
        this.k = 1;
        this.nb.setZero();
        for (int i = 0; i < numIntervals; i++) {
            processInterval(ti.getInterval(i));
            processObservation(ti.getIntervalType(i));
        }
    }

    /**
     * This method should mutate the input to account for the observation that occurred.
     *
     * @param intervalType
     *
     * @see beast.evolution.tree.birthdeath.EventType
     *
     */
    private void processObservation(EventType intervalType) {
        double lnL;

        switch (intervalType.toString()) {
            case "birth" -> {
                lnL = Math.log(birth());
                this.k += 1;
            }
            case "sample" -> {
                lnL = Math.log(psi());
                this.k -= 1;
            }
            case "occurrence" -> {
                lnL = Math.log(omega()) + this.nb.lnPGFDash1(1.0);
                double lnRp1 = Math.log(Math.exp(this.nb.getLnR()) + 1.0);
                this.nb.setLnPAndLnR(this.nb.getLnP(), lnRp1);
            }
            case "catastrophe" -> {
                int n = intervalType.getCount().getAsInt();
                double rho = p();
                lnL = (this.k - n) * Math.log(1 - rho)
                        + n * Math.log(rho)
                        + this.nb.lnPGF(1 - rho);
                this.k -= n;
                this.nb.setLnPAndLnR(Math.log(1 - rho) + this.nb.getLnP(),
                        this.nb.getLnR());

            }
            default -> throw new IllegalStateException("Unexpected value: " + intervalType.toString());
        }

        this.lnL+=lnL;
    }

    /**
     * This method should mutate the input to adjust for the interval during which there was no observation.
     *
     * @param intervalDuration the duration of time during which there was no observation
     */
    private void processInterval(double intervalDuration) {

        double lnC, lnMean, lnVariance;

        double p0Val = p0(intervalDuration); // TODO Check that this gives the correct value....
        double lnP0Dash1Val = lnP0Dash1(intervalDuration);
        double lnP0Dash2Val = lnP0Dash2(intervalDuration);
        double lnRVal = lnR(intervalDuration);
        double lnRDash1Val = lnRDash1(intervalDuration);
        double lnRDash2Val = lnRDash2(intervalDuration);

        double lnFM0, lnFM1, lnFM2;
        if (!this.nb.isZero) {
            assert this.k >= 0;
            if (this.k > 0) {

                lnFM0 = this.nb.lnPGF(p0Val) + this.k * lnRVal;
                double tmp1 = this.nb.lnPGFDash1(p0Val) + lnP0Dash1Val + this.k * lnRVal;
                double tmp2 = Math.log(k) + (k-1) * lnRVal + lnRDash1Val + this.nb.lnPGF(p0Val);
                lnFM1 = logSumExp(new double[]{tmp1, tmp2});
                tmp1 = this.nb.lnPGFDash2(p0Val) + 2 * lnP0Dash1Val + this.k * lnRVal;
                tmp2 = this.nb.lnPGFDash1(p0Val) + lnP0Dash2Val + this.k * lnRVal;
                double tmp3 = Math.log(2) + this.nb.lnPGFDash1(p0Val) + lnP0Dash1Val + Math.log(this.k) + (this.k - 1) * lnRVal + lnRDash1Val;
                double tmp4 = this.nb.lnPGF(p0Val) + Math.log(this.k) + Math.log(this.k - 1) + (this.k - 2) * lnRVal + 2 * lnRDash1Val;
                double tmp5 = this.nb.lnPGF(p0Val) + Math.log(this.k) + (this.k-1) * lnRVal + lnRDash2Val;
                lnFM2 = logSumExp(new double[]{tmp1,tmp2,tmp3,tmp4,tmp5});

            } else {

                lnFM0 = this.nb.lnPGF(p0Val);
                lnFM1 = this.nb.lnPGFDash1(p0Val) + lnP0Dash1Val;
                double tmp1 = this.nb.lnPGFDash2(p0Val) + 2 * lnP0Dash1Val;
                double tmp2 = this.nb.lnPGFDash1(p0Val) + lnP0Dash2Val;
                lnFM2 = logSumExp(new double[]{tmp1, tmp2});

            }
        } else {

            lnFM0 = lnRVal;
            lnFM1 = lnRDash1Val;
            lnFM2 = lnRDash2Val;

        }

        lnC = lnFM0;
        lnMean = lnFM1 - lnFM0;
        lnVariance = Math.log(Math.exp(lnFM2 - lnFM0) + Math.exp(lnMean) * (1 - Math.exp(lnMean)));

        this.nb.setLnMeanAndLnVariance(lnMean, lnVariance);
        this.lnL+=lnC;
    }

    private double lnR(double intervalDuration) {
        double[] tmp = odeHelpers(intervalDuration);
        double x1,x2,disc,expFact;
        x1 = tmp[0];
        x2 = tmp[1];
        disc = tmp[2];
        expFact = tmp[3];
        return Math.log(disc)
                + Math.log(expFact)
                - (2 * Math.log(birth()))
                - (2 * Math.log((x2 - 1.0) - (x1 - 1.0) * expFact));
    }

    private double lnRDash1(double intervalDuration) {
        double[] tmp = odeHelpers(intervalDuration);
        double x1,x2,disc,expFact;
        x1 = tmp[0];
        x2 = tmp[1];
        disc = tmp[2];
        expFact = tmp[3];
        return Math.log(2)
                + Math.log(1 - expFact)
                + Math.log(expFact)
                + Math.log(disc)
                - 2 * Math.log(birth())
                - 3 * Math.log((x2 - expFact * (x1 - 1.0) - 1.0));
    }

    private double lnRDash2(double intervalDuration) {
        double[] tmp = odeHelpers(intervalDuration);
        double x1,x2,disc,expFact;
        x1 = tmp[0];
        x2 = tmp[1];
        disc = tmp[2];
        expFact = tmp[3];
        return Math.log(6)
                + 2 * Math.log(1 - expFact)
                + Math.log(expFact)
                + Math.log(disc)
                - 2 * Math.log(birth())
                - 4 * Math.log((x2 - expFact * (x1 - 1.0) - 1.0));
    }

    private double p0(double intervalDuration) {
        double[] tmp = odeHelpers(intervalDuration);
        double x1 = tmp[0];
        double x2 = tmp[1];
        double expFact = tmp[3];
        return (x1 * (x2 - 1.0) - x2 * (x1 - 1.0) * expFact) / ((x2 - 1.0) - (x1 - 1.0) * expFact);
    }

    protected double p0(double intervalDuration, double z) {
        double[] tmp = odeHelpers(intervalDuration);
        double x1 = tmp[0];
        double x2 = tmp[1];
        double expFact = tmp[3];
        return (x1 * (x2 - z) - x2 * (x1 - z) * expFact) / ((x2 - z) - (x1 - z) * expFact);
    }

    protected double lnP0Dash1(double intervalDuration) {
        double[] tmp = odeHelpers(intervalDuration);
        double x1 = tmp[0];
        double x2 = tmp[1];
        double expFact = tmp[3];
        double aa = x2 - x1 * expFact;
        double bb = 1 - expFact;
        double cc = x2 * expFact - x1;
        return Math.log(cc * aa + x1 * x2 * Math.pow(bb, 2.0))
                - 2.0 * Math.log(aa - bb);
    }

    private double lnP0Dash2(double intervalDuration) {
        double[] tmp = odeHelpers(intervalDuration);
        double x1 = tmp[0];
        double x2 = tmp[1];
        double expFact = tmp[3];
        double aa = x2 - x1 * expFact;
        double bb = 1 - expFact;
        double cc = x2 * expFact - x1;
        return Math.log(2.0)
                + Math.log(bb)
                + Math.log(cc * aa + x1 * x2 * Math.pow(bb, 2.0))
                - 3.0 * Math.log(aa - bb);
    }

    private double[] odeHelpers(double intervalDuration) {
        double gamma = birth() + death() + psi() + omega();
        double discriminant = Math.pow(gamma, 2.0) - 4.0 * birth() * death();
        double sqrtDisc = Math.sqrt(discriminant);
        double x1 = (gamma - sqrtDisc) / (2 * birth());
        double x2 = (gamma + sqrtDisc) / (2 * birth());
        double expFact = Math.exp(- sqrtDisc * intervalDuration);
        return new double[]{x1, x2, discriminant, expFact};
    }

    public NegativeBinomial getNegativeBinomial() {
        return new NegativeBinomial();
    }

    public class NegativeBinomial {

        double lnMean;
        double lnVariance;
        double lnP;
        double ln1mP;
        double lnR;

        public void setZero() {
            isZero = true;
        }

        boolean isZero;

        public NegativeBinomial() {
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
            if (n < 1) {
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

        public NegativeBinomial(double mean, double variance) {
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

    /**
     * LogSumExp
     *
     * <p>This function implements a numerically safe way to take the logarithm of the sum of exponentials.</p>
     *
     * @see <a href="https://en.wikipedia.org/wiki/LogSumExp">Wikipedia page.</a>
     */
    final double logSumExp(double[] xs) {
        double xMax = Arrays.stream(xs).max().getAsDouble();
        double tmp = Arrays.stream(xs).map((x) -> Math.exp(x - xMax)).sum();
        return xMax + Math.log(tmp);
    }
}
