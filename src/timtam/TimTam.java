package timtam;


import beast.core.Citation;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeDistribution;

/**
 * Tree prior for birth-death-sampling while tracking the distribution of hidden
 * lineages. This used BirthDeathSerialSampling.java as a starting point. There
 * is a nested class, NegativeBinomial, which is the approximation used for the
 * number of unobserved lineages in the likelihood.
 *
 * @author Alexander Zarebski
 */
@Citation(value = "Zarebski AE, du Plessis L, Parag KV, Pybus OG (2022) A computationally tractable birth-death model that combines phylogenetic and epidemiological data. PLOS Computational Biology 18(2): e1009805. https://doi.org/10.1371/journal.pcbi.1009805",
        year = 2022, firstAuthorSurname = "Zarebski", DOI="10.1371/journal.pcbi.1009805")
public class TimTam extends TreeDistribution {

    final public Input<RealParameter> lambdaInput =
            new Input<>("lambda", "The birth rate. If you want to have rates that change over time you will also need the lambdaChangeTimes to be set.", (RealParameter) null);

    final public Input<RealParameter> muInput =
            new Input<>("mu", "The death rate, i.e. the rate at which individuals are removed without being observed.", (RealParameter) null);

    final public Input<RealParameter> psiInput =
            new Input<>("psi", "The sampling rate, i.e. the rate of unscheduled sequenced sampling.", (RealParameter) null);

    final public Input<RealParameter> rhoInput =
            new Input<>("rho", "The probability of sampling lineages in a scheduled sample.", (RealParameter) null);

    final public Input<RealParameter> omegaInput =
            new Input<>("omega", "The occurrence rate, i.e. the rate of unscheduled unsequenced sampling.", (RealParameter) null);

    final public Input<RealParameter> rootLengthInput =
            new Input<>("rootLength", "The length of the time between the origin and the tMRCA.", (RealParameter) null);

    final public Input<BackwardsSchedule> catastropheTimesInput = new Input<>("catastropheTimes", "the times at which a scheduled sequenced sample was attempted");

    final public Input<BackwardsPointProcess> pointsInput =
            new Input<>("points", "The times at which there was an occurrence (omega-sample) event.");

    final public Input<RealParameter> nuInput = new Input<>("nu", "the probability of unsequenced scheduled sampling", Input.Validate.OPTIONAL);

    final public Input<BackwardsSchedule> disasterTimesInput = new Input<>("disasterTimes", "the times at which a scheduled unsequenced sample was attempted", Input.Validate.OPTIONAL);

    final public Input<BackwardsCounts> disasterCountsInput = new Input<>("disasterCounts", "the size of each scheduled unsequenced sample", Input.Validate.OPTIONAL);

    final public Input<Boolean> conditionOnObservationInput = new Input<>("conditionOnObservation", "if is true then condition on sampling at least one individual (psi-sampling). The default value is true.", true);

    Tree tree;
    protected Double lambda;
    protected Double mu;
    protected Double psi;
    protected Double rho;
    protected Double omega;
    private boolean hasPositiveOmega; // because there is a decent chance someone will want to set it to zero.
    protected Double nu;
    protected Double rootLength;

    // the times at which a scheduled sequenced sample was attempted
    BackwardsSchedule catastropheTimes;

    // the times and size of any disasters
    BackwardsSchedule disasterTimes;
    BackwardsCounts disasterCounts;

    // the times at which there was an occurrence sample.
    BackwardsPointProcess points;

    boolean conditionOnObservation;

    // we use this attribute to accumulate the log-likelihood in the calculation.
    private double lnL;

    // this is used to track the number of lineages in the reconstructed tree as we process the data.
    private int k;

    // this is used to track the distribution of hidden lineages as we process the data.
    private NegativeBinomial nb = new NegativeBinomial();

    private TreeWithBackwardsPointProcess ti;
    private int numIntervals;

    // Unclear why it is necessary, but BEAST expects there to be a zero-argument
    // constructor and if there isn't one it freaks out.
    public TimTam() {

    }

    public TimTam(
            RealParameter lambda,
            RealParameter mu,
            RealParameter psi,
            RealParameter rho,
            RealParameter omega,
            RealParameter rootLength,
            BackwardsSchedule catastropheTimes,
            BackwardsPointProcess points,
            RealParameter nu,
            BackwardsSchedule disasterTimes,
            BackwardsCounts disasterCounts,
            Boolean conditionOnObservation) {
        this("timTamModel", lambda, mu, psi, rho, omega, rootLength, catastropheTimes, points, nu,
                disasterTimes, disasterCounts, conditionOnObservation);
    }

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        this.lambda = lambdaInput.get().getValue();
        this.mu = muInput.get().getValue();
        this.psi = psiInput.get().getValue();
        if (rhoInput.get() != null) {
            this.rho = rhoInput.get().getValue();
        } else {
            this.rho = null;
        }
        if (omegaInput.get() != null) {
            this.omega = omegaInput.get().getValue();
            this.hasPositiveOmega = true;
        } else {
            this.hasPositiveOmega = false;
        }
        if (nuInput.get() != null) {
            this.nu = nuInput.get().getValue();
        } else {
            this.nu = null;
        }
        this.rootLength = rootLengthInput.get().getValue();

        this.catastropheTimes = catastropheTimesInput.get();

        this.points = pointsInput.get();

        this.disasterTimes = disasterTimesInput.get();
        this.disasterCounts = disasterCountsInput.get();

        this.conditionOnObservation = conditionOnObservationInput.get();

        this.nb.setZero();

        this.tree = (Tree) treeInput.get();
        this.ti = new TreeWithBackwardsPointProcess(
                rootLengthInput.get(),
                tree,
                points,
                disasterTimes,
                disasterCounts,
                catastropheTimes);
        this.numIntervals = this.ti.getIntervalCount();
    }

    public TimTam(
            String modelName,
            RealParameter lambda,
            RealParameter mu,
            RealParameter psi,
            RealParameter rho,
            RealParameter omega,
            RealParameter rootLength,
            BackwardsSchedule catastropheTimes,
            BackwardsPointProcess points,
            RealParameter nu,
            BackwardsSchedule disasterTimes,
            BackwardsCounts disasterCounts,
            Boolean conditionOnObservation) {

        this.lambda = lambda.getValue();

        this.mu = mu.getValue();

        this.psi = psi.getValue();

        this.rho = rho.getValue();

        this.omega = omega.getValue();
        hasPositiveOmega = this.omega != null;

        this.rootLength= rootLength.getValue();

        this.catastropheTimes = catastropheTimes;

        this.points = points;

        this.nu = nu.getValue();
        if (this.nu != null) {
            nu.setBounds(0.0, 1.0);
        }
        this.disasterTimes = disasterTimes;
        this.disasterCounts = disasterCounts;

        this.conditionOnObservation = conditionOnObservation;

        this.nb.setZero();
    }

    public double birth() {
        return lambda;
    }

    public double death() {
        return mu;
    }

    public double psi() {
        return psi;
    }

    public double omega() {
        if (hasPositiveOmega) {
            return omega;
        } else {
            throw new RuntimeException("Omega was not provided hence should not be requested.");
        }
    }

    public double rho() {
        return rho;
    }

    public double nu() {
        return nu;
    }

    @Override
    public double calculateLogP() {
        calculateTreeLogLikelihood();
        return this.lnL;
    }

    /**
     * Generic likelihood calculation
     */
    public final void calculateTreeLogLikelihood() {
        // if the likelihood conditions upon the observation of the process then we need to account for this in the
        // log-likelihood.
        if (this.conditionOnObservation) {
            double probUnobserved = p0(this.ti.getTotalTimeSpan());
            this.lnL = - Math.log(1 - probUnobserved);
        } else {
            this.lnL = 0.0;
        }

        // TODO this assumes a fixed initial condition which should be moved
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
     * @param intervalType the type of observation that was made
     *
     * @see timtam.EventType
     *
     */
    private void processObservation(EventType intervalType) {
        switch (intervalType.toString()) {
            case "birth" -> {
                this.lnL += Math.log(birth());
                this.k += 1;
            }
            case "sample" -> {
                this.lnL += Math.log(psi());
                this.k -= 1;
            }
            case "occurrence" -> {
                this.lnL += Math.log(omega()) + this.nb.lnPGFDash1(1.0);
                double lnRp1 = Math.log(Math.exp(this.nb.getLnR()) + 1.0);
                this.nb.setLnPAndLnR(this.nb.getLnP(), lnRp1);
            }
            case "catastrophe" -> {
                int n = intervalType.getCount().getAsInt();
                double rho = rho();
                this.lnL += (this.k - n) * Math.log(1 - rho)
                        + n * Math.log(rho)
                        + this.nb.lnPGF(1 - rho);
                this.k -= n;
                this.nb.setLnPAndLnR(Math.log(1 - rho) + this.nb.getLnP(),
                        this.nb.getLnR());

            }
            case "disaster" -> {
                int h = intervalType.getCount().getAsInt();
                double nu = nu();
                if (h > 0) {
                    this.lnL += this.k * Math.log(1 - nu)
                            + h * Math.log(nu)
                            + this.nb.lnPGFDash(h, 1 - nu);
                } else if (h == 0) {
                    this.lnL += this.k * Math.log(1 - nu)
                            + this.nb.lnPGFDash(h, 1 - nu); // this should just be the lnPGF.
                } else {
                    throw new RuntimeException("a disaster with a negative number of cases should never happen.");
                }
                this.nb.setLnPAndLnR(
                        Math.log(1 - nu) + this.nb.getLnP(),
                        Math.log(Math.exp(this.nb.getLnR()) + h));
            }
            default -> throw new IllegalStateException("Unexpected value: " + intervalType);
        }
    }

    // this variable is just here in an attempt to resolve a memory leak...
    private double[] tmpArry = {0,0,0,0,0};

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
                tmpArry[0] = this.nb.lnPGFDash1(p0Val) + lnP0Dash1Val + this.k * lnRVal;
                tmpArry[1] = Math.log(k) + (k-1) * lnRVal + lnRDash1Val + this.nb.lnPGF(p0Val);
//                double tmp1 = this.nb.lnPGFDash1(p0Val) + lnP0Dash1Val + this.k * lnRVal;
//                double tmp2 = Math.log(k) + (k-1) * lnRVal + lnRDash1Val + this.nb.lnPGF(p0Val);
//                lnFM1 = logSumExp(new double[] {tmp1, tmp2});
                lnFM1 = logSumExp(tmpArry, 2);
                tmpArry[0] = this.nb.lnPGFDash2(p0Val) + 2 * lnP0Dash1Val + this.k * lnRVal;
                tmpArry[1] = this.nb.lnPGFDash1(p0Val) + lnP0Dash2Val + this.k * lnRVal;
                tmpArry[2] = Math.log(2) + this.nb.lnPGFDash1(p0Val) + lnP0Dash1Val + Math.log(this.k) + (this.k - 1) * lnRVal + lnRDash1Val;
                tmpArry[3] = this.nb.lnPGF(p0Val) + Math.log(this.k) + Math.log(this.k - 1) + (this.k - 2) * lnRVal + 2 * lnRDash1Val;
                tmpArry[4] = this.nb.lnPGF(p0Val) + Math.log(this.k) + (this.k-1) * lnRVal + lnRDash2Val;
//                tmp1 = this.nb.lnPGFDash2(p0Val) + 2 * lnP0Dash1Val + this.k * lnRVal;
//                tmp2 = this.nb.lnPGFDash1(p0Val) + lnP0Dash2Val + this.k * lnRVal;
//                double tmp3 = Math.log(2) + this.nb.lnPGFDash1(p0Val) + lnP0Dash1Val + Math.log(this.k) + (this.k - 1) * lnRVal + lnRDash1Val;
//                double tmp4 = this.nb.lnPGF(p0Val) + Math.log(this.k) + Math.log(this.k - 1) + (this.k - 2) * lnRVal + 2 * lnRDash1Val;
//                double tmp5 = this.nb.lnPGF(p0Val) + Math.log(this.k) + (this.k-1) * lnRVal + lnRDash2Val;
//                lnFM2 = logSumExp(new double[]{tmp1,tmp2,tmp3,tmp4,tmp5});
                lnFM2 = logSumExp(tmpArry, 5);
            } else {

                lnFM0 = this.nb.lnPGF(p0Val);
                lnFM1 = this.nb.lnPGFDash1(p0Val) + lnP0Dash1Val;
//                double tmp1 = this.nb.lnPGFDash2(p0Val) + 2 * lnP0Dash1Val;
//                double tmp2 = this.nb.lnPGFDash1(p0Val) + lnP0Dash2Val;
//                lnFM2 = logSumExp(new double[]{tmp1, tmp2});
                tmpArry[0] = this.nb.lnPGFDash2(p0Val) + 2 * lnP0Dash1Val;
                tmpArry[1] = this.nb.lnPGFDash1(p0Val) + lnP0Dash2Val;
                lnFM2 = logSumExp(tmpArry, 2);
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
        odeHelpers(intervalDuration);
//        double[] tmp = odeHelpers(intervalDuration);
//        double x1,x2,disc,expFact;
//        x1 = tmp[0];
//        x2 = tmp[1];
//        disc = tmp[2];
//        expFact = tmp[3];
//        return Math.log(disc)
//                + Math.log(expFact)
//                - (2 * Math.log(birth()))
//                - (2 * Math.log((x2 - 1.0) - (x1 - 1.0) * expFact));
        return Math.log(this.ohDiscriminant)
                + Math.log(this.ohExpFact)
                - (2 * Math.log(birth()))
                - (2 * Math.log((this.ohX2 - 1.0) - (this.ohX1 - 1.0) * this.ohExpFact));
    }

    private double lnRDash1(double intervalDuration) {
        odeHelpers(intervalDuration);
//        double[] tmp = odeHelpers(intervalDuration);
//        double x1,x2,disc,expFact;
//        x1 = tmp[0];
//        x2 = tmp[1];
//        disc = tmp[2];
//        expFact = tmp[3];
//        return Math.log(2)
//                + Math.log(1 - expFact)
//                + Math.log(expFact)
//                + Math.log(disc)
//                - 2 * Math.log(birth())
//                - 3 * Math.log((x2 - expFact * (x1 - 1.0) - 1.0));
        return Math.log(2)
                + Math.log(1 - this.ohExpFact)
                + Math.log(this.ohExpFact)
                + Math.log(this.ohDiscriminant)
                - 2 * Math.log(birth())
                - 3 * Math.log((this.ohX2 - this.ohExpFact * (this.ohX1 - 1.0) - 1.0));
    }

    private double lnRDash2(double intervalDuration) {
        odeHelpers(intervalDuration);
//        double[] tmp = odeHelpers(intervalDuration);
//        double x1,x2,disc,expFact;
//        x1 = tmp[0];
//        x2 = tmp[1];
//        disc = tmp[2];
//        expFact = tmp[3];
//        return Math.log(6)
//                + 2 * Math.log(1 - expFact)
//                + Math.log(expFact)
//                + Math.log(disc)
//                - 2 * Math.log(birth())
//                - 4 * Math.log((x2 - expFact * (x1 - 1.0) - 1.0));
        return Math.log(6)
                + 2 * Math.log(1 - this.ohExpFact)
                + Math.log(this.ohExpFact)
                + Math.log(this.ohDiscriminant)
                - 2 * Math.log(birth())
                - 4 * Math.log((this.ohX2 - this.ohExpFact * (this.ohX1 - 1.0) - 1.0));
    }

    double p0(double intervalDuration) {
        odeHelpers(intervalDuration);
        return (this.ohX1 * (this.ohX2 - 1.0) - this.ohX2 * (this.ohX1 - 1.0) * this.ohExpFact) / ((this.ohX2 - 1.0) - (this.ohX1 - 1.0) * this.ohExpFact);
//        double[] tmp = odeHelpers(intervalDuration);
//        double x1 = tmp[0];
//        double x2 = tmp[1];
//        double expFact = tmp[3];
//        return (x1 * (x2 - 1.0) - x2 * (x1 - 1.0) * expFact) / ((x2 - 1.0) - (x1 - 1.0) * expFact);
    }

    protected double p0(double intervalDuration, double z) {
        odeHelpers(intervalDuration);
        return (this.ohX1 * (this.ohX2 - z) - this.ohX2 * (this.ohX1 - z) * this.ohExpFact) / ((this.ohX2 - z) - (this.ohX1 - z) * this.ohExpFact);
//        double[] tmp = odeHelpers(intervalDuration);
//        double x1 = tmp[0];
//        double x2 = tmp[1];
//        double expFact = tmp[3];
//        return (x1 * (x2 - z) - x2 * (x1 - z) * expFact) / ((x2 - z) - (x1 - z) * expFact);
    }

    protected double lnP0Dash1(double intervalDuration) {
        odeHelpers(intervalDuration);
        double aa = this.ohX2 - this.ohX1 * this.ohExpFact;
        double bb = 1 - this.ohExpFact;
        double cc = this.ohX2 * this.ohExpFact - this.ohX1;
        return Math.log(cc * aa + this.ohX1 * this.ohX2 * Math.pow(bb, 2.0))
                - 2.0 * Math.log(aa - bb);
//        double[] tmp = odeHelpers(intervalDuration);
//        double x1 = tmp[0];
//        double x2 = tmp[1];
//        double expFact = tmp[3];
//        double aa = x2 - x1 * expFact;
//        double bb = 1 - expFact;
//        double cc = x2 * expFact - x1;
//        return Math.log(cc * aa + x1 * x2 * Math.pow(bb, 2.0))
//                - 2.0 * Math.log(aa - bb);
    }

    private double lnP0Dash2(double intervalDuration) {
        odeHelpers(intervalDuration);
        double aa = this.ohX2 - this.ohX1 * this.ohExpFact;
        double bb = 1 - this.ohExpFact;
        double cc = this.ohX2 * this.ohExpFact - this.ohX1;
        return Math.log(2.0)
                + Math.log(bb)
                + Math.log(cc * aa + this.ohX1 * this.ohX2 * Math.pow(bb, 2.0))
                - 3.0 * Math.log(aa - bb);
        // double[] tmp = odeHelpers(intervalDuration);
        // double x1 = tmp[0];
        // double x2 = tmp[1];
        // double expFact = tmp[3];
        // double aa = x2 - x1 * expFact;
        // double bb = 1 - expFact;
        // double cc = x2 * expFact - x1;
        // return Math.log(2.0)
        //         + Math.log(bb)
        //         + Math.log(cc * aa + x1 * x2 * Math.pow(bb, 2.0))
        //         - 3.0 * Math.log(aa - bb);
    }

    private double ohX1, ohX2, ohDiscriminant, ohExpFact;

    private void odeHelpers(double intervalDuration) {
        double gamma;
        if (hasPositiveOmega) {
            gamma = birth() + death() + psi() + omega();
        } else {
            gamma = birth() + death() + psi();
        }
        this.ohDiscriminant = Math.pow(gamma, 2.0) - 4.0 * birth() * death();
        double sqrtDisc = Math.sqrt(this.ohDiscriminant);
        this.ohX1 = (gamma - sqrtDisc) / (2 * birth());
        this.ohX2 = (gamma + sqrtDisc) / (2 * birth());
        this.ohExpFact = Math.exp(- sqrtDisc * intervalDuration);
    }

    public NegativeBinomial getNegativeBinomial() {
        return this.nb;
    }

    public NegativeBinomial getNewNegativeBinomial() {
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
    final static double logSumExp(double[] xs) {
        return logSumExp(xs, xs.length);
    }

    /**
     * This is a nitty gritty log-sum-exp which is useful when fine tuning memory usage by avoiding streams.
     *
     * @param xs array
     * @param n number of leading entries to use
     */
    final static double logSumExp(double[] xs, int n) {
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
       // double xMax = Arrays.stream(xs).limit(n).max().getAsDouble();
       // double tmp = Arrays.stream(xs).limit(n).map((x) -> Math.exp(x - xMax)).sum();
       // return xMax + Math.log(tmp);
    }

    @Override
    protected boolean requiresRecalculation() {
        // We need to tell beast that the likelihood needs to be recalculated each time one of the inputs corresponding
        // to a parameter changes. Since the schedule and points are not estimated we do not care about them here.
        return true;
//        return super.requiresRecalculation()
//                || lambdaInput.get().somethingIsDirty()
//                || muInput.get().somethingIsDirty()
//                || psiInput.get().somethingIsDirty()
//                || rhoInput.get().somethingIsDirty()
//                || omegaInput.get().somethingIsDirty()
//                || nuInput.get().somethingIsDirty()
//                || rootLengthInput.get().somethingIsDirty();
    }
}
