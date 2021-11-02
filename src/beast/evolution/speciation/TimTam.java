package beast.evolution.speciation;


import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.birthdeath.EventType;
import beast.evolution.tree.birthdeath.TreeWithPointProcess;
import beast.evolution.tree.coalescent.TreeIntervals;
import cern.jet.random.NegativeBinomial;

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

    boolean hasFinalSample = false;

    // length of the edge from origin to MRCA node.
    RealParameter rootLength;

    // the times at which there was an occurrence sample.
    TraitSet points;

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
            boolean hasFinalSample,
            RealParameter rootLength,
            TraitSet points) {
        //Type units) {

        this("timTamModel", lambda, mu, psi, p, omega, hasFinalSample, rootLength, points);
    }

    final public Input<RealParameter> lambdaInput = new Input<>("lambda", "the birth rate of new infections");
    final public Input<RealParameter> muInput = new Input<>("mu", "the death rate");
    final public Input<RealParameter> psiInput = new Input<>("psi", "the sampling rate");
    final public Input<RealParameter> pInput = new Input<>("p", "the probability of sampling extant lineages");
    final public Input<RealParameter> omegaInput = new Input<>("omega", "the occurrence rate");
    final public Input<Boolean> hasFinalSampleInput = new Input<>("hasFinalSample", "boolean for if there was a scheduled sample at the present");
    final public Input<RealParameter> rootLengthInput = new Input<>("rootLength", "the length of the edge between the origin and the MRCA");
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

        this.hasFinalSample = hasFinalSampleInput.get();

        this.rootLength= rootLengthInput.get();
        rootLength.setBounds(0.0, Double.POSITIVE_INFINITY);

        this.points = pointsInput.get();
    }

    public TimTam(
            String modelName,
            RealParameter lambda,
            RealParameter mu,
            RealParameter psi,
            RealParameter p,
            RealParameter omega,
            boolean hasFinalSample,
            RealParameter rootLength,
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

        this.hasFinalSample = hasFinalSample;

        this.rootLength= rootLength;
        rootLength.setBounds(0.0, Double.POSITIVE_INFINITY);

        this.points = points;
    }

    /**
     * @param b   birth rate
     * @param d   death rate
     * @param p   proportion sampled at final time point
     * @param psi rate of sampling per lineage per unit time
     * @param t   time
     * @return the probability of no sampled descendants after time, t
     */
    public static double p0(double b, double d, double p, double psi, double t) {
        double c1 = c1(b, d, psi);
        double c2 = c2(b, d, p, psi);

        double expc1trc2 = Math.exp(-c1 * t) * (1.0 - c2);

        return (b + d + psi + c1 * ((expc1trc2 - (1.0 + c2)) / (expc1trc2 + (1.0 + c2)))) / (2.0 * b);
    }

    public static double q(double b, double d, double p, double psi, double t) {
        double c1 = c1(b, d, psi);
        double c2 = c2(b, d, p, psi);
//        double res = 2.0 * (1.0 - c2 * c2) + Math.exp(-c1 * t) * (1.0 - c2) * (1.0 - c2) + Math.exp(c1 * t) * (1.0 + c2) * (1.0 + c2);
        double res = c1 * t + 2.0 * Math.log(Math.exp(-c1 * t) * (1.0 - c2) + (1.0 + c2)); // operate directly in logspace, c1 * t too big
        return res;
    }

    private static double c1(double b, double d, double psi) {
        return Math.abs(Math.sqrt(Math.pow(b - d - psi, 2.0) + 4.0 * b * psi));
    }

    private static double c2(double b, double d, double p, double psi) {
        return -(b - d - 2.0 * b * p - psi) / c1(b, d, psi);
    }


    public double p0(double t) {
        return p0(birth(), death(), p(), psi(), t);
    }

    public double q(double t) {
        return q(birth(), death(), p(), psi(), t);
    }

    private double c1() {
        return c1(birth(), death(), psi());
    }

    private double c2() {
        return c2(birth(), death(), p(), psi());
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

    /**
     * @return the proportion of population sampled at final sample, or zero if there is no final sample
     */
    public double p() {

//        if (mask != null) return mask.p.getValue(0);
        return hasFinalSample ? p.getValue(0) : 0;
    }

    // The mask does not affect the following two methods

    public double x0() {
        return rootLength.getValue(0);
    }

    @Override
    public double calculateLogP() {
        logP = calculateTreeLogLikelihood((Tree) treeInput.get());
        return logP;
    }

    /**
     * Generic likelihood calculation
     *
     * @param tree the tree to calculate likelihood of
     * @return log-likelihood of density
     */
    public final double calculateTreeLogLikelihood(Tree tree) {
        System.out.println("the calculateTreeLogLikelihood method has been called...");

        TreeWithPointProcess ti = new TreeWithPointProcess(rootLength, tree, points);
        int numIntervals = ti.getIntervalCount();

        double lnP = 0.0;

        // TODO this assumes a fixed initial condition which should be movied
        //  into something that gets specified by the client.
        int k = 1;
        NegativeBinomial hiddenLineages = new NegativeBinomial();
        hiddenLineages.setZero(true);

        for (int i = 0; i <= numIntervals; i++) {
            processInterval(ti.getInterval(i), k, hiddenLineages, lnP);
            processObservation(ti.getIntervalType(i), k, hiddenLineages, lnP);
        }

        return lnP;

        // it is impossible for the origin to be closer to the present than the
        // depth of the tree so if this is the case we can return negative
        // infinity for the log-likelihood immediately.
//        if (x0() < tree.getRoot().getHeight()) {
//            return Double.NEGATIVE_INFINITY;
//        }

        // extant leaves
//        int n = 0;
        // extinct leaves
//        int m = 0;

//        for (int i = 0; i < tree.getLeafNodeCount(); i++) {
//            Node node = tree.getNode(i);
//            if (node.getHeight() == 0.0) {
//                n += 1;
//            } else {
//                m += 1;
//            }
//        }
//
//        if (!hasFinalSample && n < 1) {
//            throw new RuntimeException(
//                    "For sampling-through-time model there must be at least one tip at time zero.");
//        }
//
//        double b = birth();
//        double p = p();
//
//        double logL;
//        logL = -q(x0());
//        if (hasFinalSample) {
//            logL += n * Math.log(4.0 * p);
//        }
//        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
//            double x = tree.getNode(tree.getLeafNodeCount() + i).getHeight();
//            logL += Math.log(b) - q(x);
//
//            //System.out.println("internalNodeLogL=" + Math.log(b / q(x)));
//
//        }
//        for (int i = 0; i < tree.getLeafNodeCount(); i++) {
//            double y = tree.getNode(i).getHeight();
//
//            if (y > 0.0) {
//                logL += Math.log(psi()) + q(y);
//
//                //System.out.println("externalNodeLogL=" + Math.log(psi() * (r() + (1.0 - r()) * p0(y)) * q(y)));
//
//            } else if (!hasFinalSample) {
//                //handle condition ending on final tip in sampling-through-time-only situation
//                logL += Math.log(psi()) + q(y);
////                System.out.println("externalNodeLogL=" + Math.log(psi() * q(y)));
//
//            }
//        }
//
//        return logL;
    }

    /**
     * This method should mutate the input to account for the observation that occurred.
     *
     * @param intervalType
     * @param k
     * @param hiddenLineages
     * @param lnP
     */
    private void processObservation(EventType intervalType, int k, NegativeBinomial hiddenLineages, double lnP) {
    }

    /**
     * This method should mutate the input
     * @param interval
     * @param k
     * @param hiddenLineages
     * @param lnP
     */
    private void processInterval(double interval, int k, NegativeBinomial hiddenLineages, double lnP) {
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

        public void setZero(boolean zero) {
            isZero = zero;
        }

        boolean isZero;

        public NegativeBinomial() {
        }

        /**
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
        }

        public double getLnP() {
            return lnP;
        }

        public double getLnR() {
            return lnR;
        }

        public void setLnPAndLnR(double lnP, double lnR) {
            if (!this.isZero) {
                this.lnP = lnP;
                this.ln1mP = Math.log(1 - Math.exp(lnP));
                this.lnR = lnR;
                updateLnMeanAndLnVariance();
            } else {
                throw new RuntimeException("r and p of negative binomial not defined for zero distribution.");
            }
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
}
