package timtam;


import beast.core.Citation;
import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeDistribution;

import java.util.Arrays;
import java.util.Objects;
import java.util.OptionalInt;

/**
 * Tree prior for birth-death-sampling while tracking the distribution of hidden
 * lineages. This used BirthDeathSerialSampling.java as a starting point. There
 * is a nested class, NegativeBinomial, which is the approximation used for the
 * number of unobserved lineages in the likelihood.
 *
 * Time is measured backwards from the last sequenced sample (which is treated
 * as having been collected at time zero).
 *
 * @author Alexander Zarebski
 */
@Citation(value = "Zarebski AE, du Plessis L, Parag KV, Pybus OG (2022) A computationally tractable birth-death model that combines phylogenetic and epidemiological data. PLOS Computational Biology 18(2): e1009805. https://doi.org/10.1371/journal.pcbi.1009805",
        year = 2022, firstAuthorSurname = "Zarebski", DOI="10.1371/journal.pcbi.1009805")
public class TimTam extends TreeDistribution {

    public Input<RealParameter> lambdaInput =
            new Input<>("lambda", "The birth rate. If you want to have rates that change over time you will also need the lambdaChangeTimes to be set.", (RealParameter) null);

    public Input<RealParameter> lambdaChangeTimesInput =
            new Input<>("lambdaChangeTimes", "The times at which the value of lambda changes. These should be given as backwards times treating the final observation *in the tree* as the present (time zero). If lambda is constant then this parameter can safely be left as the default null value.", (RealParameter) null);

    public Input<RealParameter> muInput =
            new Input<>("mu", "The death rate, i.e. the rate at which individuals are removed without being observed.", (RealParameter) null);

    public Input<RealParameter> psiInput =
            new Input<>("psi", "The sampling rate, i.e. the rate of unscheduled sequenced sampling.", (RealParameter) null);

    public Input<RealParameter> rhoInput =
            new Input<>("rho", "The probability of sampling lineages in a scheduled sample.", (RealParameter) null);

    public Input<RealParameter> omegaInput =
            new Input<>("omega", "The occurrence rate, i.e. the rate of unscheduled unsequenced sampling. Default value of zero (no occurrence sampling).", (RealParameter) null);

    public Input<RealParameter> originTimeInput =
            new Input<>("originTime",
                    "The (backwards) time of the origin (relative to the most recent observation in the tree which is considered as having occurred at time zero).",
                    (RealParameter) null);

    public Input<RealParameter> catastropheTimesInput = new Input<>("catastropheTimes", "the times at which a scheduled sequenced sample was attempted", Input.Validate.OPTIONAL);

    public Input<RealParameter> occurrenceTimesInput =
            new Input<>("occurrenceTimes",
                    "The times at which there was an occurrence (omega-sample) event. This should be entered as backwards time relative to the final time in the tree which is treated as the present (time zero). The default variable for this is null meaning no occurrence observations, leaving this empty does not mean that the omega rate is assumed to be zero.");

    public Input<RealParameter> nuInput =
            new Input<>("nu", "the probability of unsequenced scheduled sampling", Input.Validate.OPTIONAL);

    public Input<RealParameter> disasterTimesInput =
            new Input<>("disasterTimes",
                    "The times at which a scheduled unsequenced sample was attempted. This should be entered as backwards time relative to the final time in the tree which is treated as the present (time zero). The default variable for this is null indicating that no scheduled unsequenced samples were attempted.",
                    (RealParameter) null);

    public Input<IntegerParameter> disasterSizesInput =
            new Input<>("disasterSizes",
                    "the size of each scheduled unsequenced sample",
                    (IntegerParameter) null);

    public Input<Boolean> conditionOnObservationInput = new Input<>("conditionOnObservation", "if is true then condition on sampling at least one individual (psi-sampling). The default value is true.", true);

    // we specify a threshold below which two times are considered equal.
    private final double timeEpsilon = 0.00001;

    Tree tree;
    protected Double[] lambda;
    protected Double[] lambdaChangeTimes;
    protected Double mu;
    protected Double psi;
    protected Double rho;
    protected Double omega;
    protected Double nu;
    protected Double originTime;

    private Double[] rateChangeTimes;
    protected int numRateChanges;

    protected Double[] catastropheTimes;
    protected Integer[] catastropheSizes;
    private Integer totalCatastropheSizes;

    // A disaster is a scheduled sample of lineages where sampled lineages are *not* sequenced so do not appear in the
    // reconstructed tree. Typically, these data will form a time series.
    protected Double[] disasterTimes;
    protected Integer[] disasterSizes;

    // the times at which there was an occurrence sample measured in backwards time with the final tip in the tree being
    // used as zero. An occurrence is an unscheduled and unsequenced sample. These data form a point process of events.
    protected Double[] occurrenceTimes;

    boolean conditionOnObservation;

    // we use this attribute to accumulate the log-likelihood in the calculation.
    private double lnL;

    // this is used to track the number of lineages in the reconstructed tree as we process the data.
    private int k;

    // this is used to track the distribution of hidden lineages as we process the data.
    private final TimTamNegBinom nb = new TimTamNegBinom();



    // this is the time from the origin until the time of the final observation.
    private double timeFromOriginToFinalDatum;

    // this is the number of intervals of time that need to be considered when evaluating the likelihood. An interval
    // could be due to a change in a rate parameter or an observation, e.g. a birth or an occurrence or a disaster.
    private int numTimeIntervals;

    // this array holds the reason each interval terminated, this could be because there was an observation or a rate
    // change.
    private TimTamIntervalTerminator[] intervalTerminators;

    // this array holds the length of each interval.
    private double[] intervalDuration;

    // Unclear why it is necessary, but BEAST expects there to be a zero-argument
    // constructor and if there isn't one it freaks out.
    public TimTam() {

    }

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        this.tree = (Tree) treeInput.get();

        updateRateAndProbParams();

        this.originTime = originTimeInput.get().getValue();

        if (catastropheTimesInput.get() != null) {
            this.catastropheTimes = catastropheTimesInput.get().getValues();
            this.catastropheSizes = new Integer[this.catastropheTimes.length];
            measureCatastrophes();
        } else {
            this.catastropheTimes = new Double[] {};
            this.catastropheSizes = new Integer[] {};
            this.totalCatastropheSizes = 0;
        }

        this.occurrenceTimes = occurrenceTimesInput.get() != null ? occurrenceTimesInput.get().getValues() : new Double[] {};

        if (disasterTimesInput.get() != null) {
            this.disasterTimes = disasterTimesInput.get().getValues();
            this.disasterSizes = disasterSizesInput.get().getValues();
        } else {
            this.disasterTimes = new Double[] {};
            this.disasterSizes = new Integer[] {};
        }

        this.conditionOnObservation = conditionOnObservationInput.get();

        this.nb.setZero();

        // We need to be careful in counting the number of time intervals because there may be multiple
        // tree nodes counted as a single catastrophe.
        this.numTimeIntervals =
                numRateChanges +
                this.disasterTimes.length + this.occurrenceTimes.length +
                this.tree.getNodesAsArray().length - this.totalCatastropheSizes + this.catastropheTimes.length;
        this.intervalTerminators = new TimTamIntervalTerminator[this.numTimeIntervals];
        int iTx = 0;

        for (Double rateChangeTime : rateChangeTimes) {
            this.intervalTerminators[iTx] =
                    new TimTamIntervalTerminator("rateChange", rateChangeTime, OptionalInt.empty());
            iTx++;
        }

        for (Double occurrenceTime : this.occurrenceTimes) {
            this.intervalTerminators[iTx] =
                    new TimTamIntervalTerminator("occurrence", occurrenceTime, OptionalInt.empty());
            iTx++;
        }

        for (int ix = 0; ix < this.catastropheTimes.length; ix++) {
            this.intervalTerminators[iTx] =
                    new TimTamIntervalTerminator("catastrophe", this.catastropheTimes[ix], OptionalInt.of(this.catastropheSizes[ix]));
            iTx++;
        }

        for (int ix = 0; ix < this.disasterTimes.length; ix++) {
            this.intervalTerminators[iTx] =
                    new TimTamIntervalTerminator("disaster", this.disasterTimes[ix], OptionalInt.of(this.disasterSizes[ix]));
            iTx++;
        }

        for (Node node : this.tree.getNodesAsArray()) {
            if (isUnscheduledTreeNode(node)) {
                this.intervalTerminators[iTx] =
                        new TimTamIntervalTerminator(node.isLeaf() ? "sample" : "birth", node.getHeight(), OptionalInt.empty());
                iTx++;
            }
        }
        Arrays.sort(this.intervalTerminators);

        this.intervalDuration = new double[this.numTimeIntervals];
        this.intervalDuration[0] = originTime - this.intervalTerminators[0].getBwdTime();
        for (int ix = 1; ix < this.intervalTerminators.length; ix++) {
            this.intervalDuration[ix] = this.intervalTerminators[ix-1].getBwdTime() - this.intervalTerminators[ix].getBwdTime();
        }

        this.timeFromOriginToFinalDatum = 0;
        for (Double intDur : intervalDuration) {
            this.timeFromOriginToFinalDatum += intDur;
        }
    }

    /**
     * Predicate for being an unscheduled node.
     */
    private boolean isUnscheduledTreeNode(Node node) {
        if (!node.isLeaf()) {
            return true;
        }
        for (Double ct : this.catastropheTimes) {
            if (Math.abs(node.getHeight() - ct) < this.timeEpsilon) {
                return false;
            }
        }
        return true;
    }

    /**
     * Loop over the tree and count the number of tips in each catastrophe. The member variables of the tree and the
     * timing of catastrophes will allow this to be done.
     */
    private void measureCatastrophes() {

        Node[] nodes = this.tree.getNodesAsArray();
        Double cT;
        int cS;

        this.totalCatastropheSizes = 0;
        for (int ix = 0; ix < this.catastropheTimes.length; ix++) {
            cT = this.catastropheTimes[ix];
            cS = 0;
            for (int jx = 0; jx < this.tree.getLeafNodeCount(); jx++) {
                if (Math.abs(nodes[jx].getHeight() - cT) < this.timeEpsilon) {
                    cS++;
                }
            }
            this.catastropheSizes[ix] = cS;
            this.totalCatastropheSizes += cS;
        }
    }

    public double birth(double bwdTime) {
        if (this.lambda.length > 1) {
            double cBwdTime = Double.POSITIVE_INFINITY;
            double cLambda = this.lambda[0];
            double nBwdTime = this.lambdaChangeTimes[0];
            int ix = 0;
            while (true) {
                if (cBwdTime > bwdTime & bwdTime > nBwdTime) {
                    return cLambda;
                } else {
                    ix+=1;
                    cBwdTime = this.lambdaChangeTimes[ix-1];
                    cLambda = this.lambda[ix];
                    nBwdTime = this.lambdaChangeTimes[ix];
                }
            }
        } else {
            return this.lambda[0];
        }
    }

    public double death() {
        return mu;
    }

    public double psi() {
        return psi;
    }

    public double omega() {
        return this.omega;
    }

    public double rho() {
        return rho;
    }

    public double nu() {
        return nu;
    }

    @Override
    public double calculateLogP() {
        updateRateAndProbParams();
        calculateTreeLogLikelihood();
        return this.lnL;
    }

    /**
     * Update the (primitive) local parameters to match the input.
     */
    private void updateRateAndProbParams() {
        this.lambda = lambdaInput.get().getValues();
        this.lambdaChangeTimes =
                (lambdaChangeTimesInput.get() != null) ? lambdaChangeTimesInput.get().getValues() : new Double[]{};
        this.mu = muInput.get().getValue();
        this.psi = psiInput.get().getValue();
        this.rho = (rhoInput.get() != null) ? rhoInput.get().getValue() : null;
        this.omega = (omegaInput.get() != null) ? omegaInput.get().getValue() : 0.0;
        this.nu = (nuInput.get() != null) ? nuInput.get().getValue() : null;

        this.rateChangeTimes = this.lambdaChangeTimes;
        this.numRateChanges = this.rateChangeTimes.length;
    }

    /**
     * Generic likelihood calculation
     */
    public final void calculateTreeLogLikelihood() {
        // if the likelihood conditions upon the observation of the process then we need to account for this in the
        // log-likelihood.
        if (this.conditionOnObservation) {
            double probUnobserved = p0(this.timeFromOriginToFinalDatum, 0.0, 1.0);
            this.lnL = - Math.log(1 - probUnobserved);
        } else {
            this.lnL = 0.0;
        }

        // This assumes a fixed initial condition (i.e. there is a single phylogeny involved). It should really be
        // moved into an input so that multiple introductions can be handled elegantly.
        this.k = 1;

        this.nb.setZero();

        for (int ix = 0; ix < this.numTimeIntervals; ix++) {
            processInterval(this.intervalDuration[ix], this.intervalTerminators[ix].getBwdTime());
            processObservation(this.intervalTerminators[ix]);
        }
    }

    /**
     * This method should mutate the attributes of this object to account for
     * the observation that occurred.
     *
     * @implNote This method is implemented using if-else rather than a switch
     * because the switch syntax has changed between Java versions so this
     * guards against changes in versions.
     *
     * @param intTerminator the reason the interval of time ended.
     *
     * @see TimTamIntervalTerminator
     *
     */
    private void processObservation(TimTamIntervalTerminator intTerminator) {
        String intTypeStr = intTerminator.getType();
        Double obsBwdTime = intTerminator.getBwdTime();

        if (Objects.equals(intTypeStr, "birth")) {
            this.lnL += Math.log(birth(obsBwdTime));
            this.k += 1;
        } else if (Objects.equals(intTypeStr, "sample")) {
            this.lnL += Math.log(psi());
            this.k -= 1;
        } else if (Objects.equals(intTypeStr, "occurrence")) {
            this.lnL += Math.log(omega()) + this.nb.lnPGFDash1(1.0);
            double lnRp1 = Math.log(Math.exp(this.nb.getLnR()) + 1.0);
            this.nb.setLnPAndLnR(this.nb.getLnP(), lnRp1);
        } else if (Objects.equals(intTypeStr, "catastrophe")) {
            int n = intTerminator.getCount();
            double rho = rho();
            this.lnL += (this.k - n) * Math.log(1 - rho)
                    + n * Math.log(rho)
                    + this.nb.lnPGF(1 - rho);
            this.k -= n;
            this.nb.setLnPAndLnR(Math.log(1 - rho) + this.nb.getLnP(),
                    this.nb.getLnR());
        } else if (Objects.equals(intTypeStr, "disaster")) {
            int h = intTerminator.getCount();
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
        } else {
            throw new IllegalStateException("Unexpected value: " + intTerminator.getType() + "\n\tPlease look at the TimTamIntervalTerminator class to see the type of intervals that TimTam recognises.");
        }
    }

    // this variable is just here in an attempt to resolve a memory leak...
    private final double[] tmpArry = {0,0,0,0,0};

    /**
     * This method should mutate the input to adjust for the interval during which there was no observation.
     *
     * @param intervalDuration the duration of time during which there was no observation
     */
    private void processInterval(double intervalDuration, double bwdTimeIntervalEnd) {

        double lnC, lnMean, lnVariance;

        double p0Val = p0(intervalDuration, bwdTimeIntervalEnd, 1.0);
        double lnP0Dash1Val = lnP0Dash1(intervalDuration, bwdTimeIntervalEnd);
        double lnP0Dash2Val = lnP0Dash2(intervalDuration, bwdTimeIntervalEnd);
        double lnRVal = lnR(intervalDuration, bwdTimeIntervalEnd);
        double lnRDash1Val = lnRDash1(intervalDuration, bwdTimeIntervalEnd);
        double lnRDash2Val = lnRDash2(intervalDuration, bwdTimeIntervalEnd);

        double lnFM0, lnFM1, lnFM2;
        if (!this.nb.isZero) {
            assert this.k >= 0;
            if (this.k > 0) {

                lnFM0 = this.nb.lnPGF(p0Val) + this.k * lnRVal;
                tmpArry[0] = this.nb.lnPGFDash1(p0Val) + lnP0Dash1Val + this.k * lnRVal;
                tmpArry[1] = Math.log(k) + (k-1) * lnRVal + lnRDash1Val + this.nb.lnPGF(p0Val);
                lnFM1 = logSumExp(tmpArry, 2);
                tmpArry[0] = this.nb.lnPGFDash2(p0Val) + 2 * lnP0Dash1Val + this.k * lnRVal;
                tmpArry[1] = this.nb.lnPGFDash1(p0Val) + lnP0Dash2Val + this.k * lnRVal;
                tmpArry[2] = Math.log(2) + this.nb.lnPGFDash1(p0Val) + lnP0Dash1Val + Math.log(this.k) + (this.k - 1) * lnRVal + lnRDash1Val;
                tmpArry[3] = this.nb.lnPGF(p0Val) + Math.log(this.k) + Math.log(this.k - 1) + (this.k - 2) * lnRVal + 2 * lnRDash1Val;
                tmpArry[4] = this.nb.lnPGF(p0Val) + Math.log(this.k) + (this.k-1) * lnRVal + lnRDash2Val;
                lnFM2 = logSumExp(tmpArry, 5);
            } else {

                lnFM0 = this.nb.lnPGF(p0Val);
                lnFM1 = this.nb.lnPGFDash1(p0Val) + lnP0Dash1Val;
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

    private double lnR(double intervalDuration, double bwdTimeIntervalEnd) {
        odeHelpers(intervalDuration, bwdTimeIntervalEnd);
        return Math.log(this.ohDiscriminant)
                + Math.log(this.ohExpFact)
                - (2 * Math.log(birth(bwdTimeIntervalEnd)))
                - (2 * Math.log((this.ohX2 - 1.0) - (this.ohX1 - 1.0) * this.ohExpFact));
    }

    private double lnRDash1(double intervalDuration, double bwdTimeIntervalEnd) {
        odeHelpers(intervalDuration, bwdTimeIntervalEnd);
        return Math.log(2)
                + Math.log(1 - this.ohExpFact)
                + Math.log(this.ohExpFact)
                + Math.log(this.ohDiscriminant)
                - 2 * Math.log(birth(bwdTimeIntervalEnd))
                - 3 * Math.log((this.ohX2 - this.ohExpFact * (this.ohX1 - 1.0) - 1.0));
    }

    private double lnRDash2(double intervalDuration, double bwdTimeIntervalEnd) {
        odeHelpers(intervalDuration, bwdTimeIntervalEnd);
        return Math.log(6)
                + 2 * Math.log(1 - this.ohExpFact)
                + Math.log(this.ohExpFact)
                + Math.log(this.ohDiscriminant)
                - 2 * Math.log(birth(bwdTimeIntervalEnd))
                - 4 * Math.log((this.ohX2 - this.ohExpFact * (this.ohX1 - 1.0) - 1.0));
    }

    protected double p0(double intervalDuration, double bwdTimeIntervalEnd, double z) {
        odeHelpers(intervalDuration, bwdTimeIntervalEnd);
        return (this.ohX1 * (this.ohX2 - z) - this.ohX2 * (this.ohX1 - z) * this.ohExpFact) / ((this.ohX2 - z) - (this.ohX1 - z) * this.ohExpFact);
    }

    protected double lnP0Dash1(double intervalDuration, double bwdTimeIntervalEnd) {
        odeHelpers(intervalDuration, bwdTimeIntervalEnd);
        double aa = this.ohX2 - this.ohX1 * this.ohExpFact;
        double bb = 1 - this.ohExpFact;
        double cc = this.ohX2 * this.ohExpFact - this.ohX1;
        return Math.log(cc * aa + this.ohX1 * this.ohX2 * Math.pow(bb, 2.0))
                - 2.0 * Math.log(aa - bb);
    }

    private double lnP0Dash2(double intervalDuration, double bwdTimeIntervalEnd) {
        odeHelpers(intervalDuration, bwdTimeIntervalEnd);
        double aa = this.ohX2 - this.ohX1 * this.ohExpFact;
        double bb = 1 - this.ohExpFact;
        double cc = this.ohX2 * this.ohExpFact - this.ohX1;
        return Math.log(2.0)
                + Math.log(bb)
                + Math.log(cc * aa + this.ohX1 * this.ohX2 * Math.pow(bb, 2.0))
                - 3.0 * Math.log(aa - bb);
    }

    private double ohX1, ohX2, ohDiscriminant, ohExpFact;

    private void odeHelpers(double intervalDuration, double bwdTimeIntervalEnd) {
        double gamma = birth(bwdTimeIntervalEnd) + death() + psi() + omega();
        this.ohDiscriminant = Math.pow(gamma, 2.0) - 4.0 * birth(bwdTimeIntervalEnd) * death();
        double sqrtDisc = Math.sqrt(this.ohDiscriminant);
        this.ohX1 = (gamma - sqrtDisc) / (2 * birth(bwdTimeIntervalEnd));
        this.ohX2 = (gamma + sqrtDisc) / (2 * birth(bwdTimeIntervalEnd));
        this.ohExpFact = Math.exp(- sqrtDisc * intervalDuration);
    }

    public TimTamNegBinom getTimTamNegBinom() {
        return this.nb;
    }

    /**
     * LogSumExp
     *
     * <p>This function implements a numerically safe way to take the logarithm of the sum of exponentials.</p>
     *
     * @see <a href="https://en.wikipedia.org/wiki/LogSumExp">Wikipedia page.</a>
     */
    final double logSumExp(double[] xs) {
        return logSumExp(xs, xs.length);
    }

    /**
     * This is a nitty gritty log-sum-exp which is useful when fine tuning memory usage by avoiding streams.
     *
     * @param xs array
     * @param n number of leading entries to use
     */
    final double logSumExp(double[] xs, int n) {
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

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }
}
