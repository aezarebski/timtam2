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
            new Input<>("lambda", "The birth rate, i.e. the rate at which an infected individual spawns another infected individual. If you want to have rates that change over time you will also need the lambdaChangeTimes to be set.", (RealParameter) null);

    public Input<RealParameter> lambdaChangeTimesInput =
            new Input<>("lambdaChangeTimes", "The times at which the value of lambda changes. These should be given as backwards times treating the final observation *in the tree* as the present (time zero). If lambda is constant then this parameter can safely be left as the default null value.", (RealParameter) null);

    public Input<RealParameter> muInput =
            new Input<>("mu", "The death rate, i.e. the rate at which individuals are removed without being observed. If you want to have rates that change over time you will also need the muChangeTimes to be set.", (RealParameter) null);

    public Input<RealParameter> muChangeTimesInput =
            new Input<>("muChangeTimes", "The times at which the value of mu changes. These should be given as backwards times treating the final observation *in the tree* as the present (time zero). If mu is constant then this parameter can safely be left as the default null value.", (RealParameter) null);

    public Input<RealParameter> psiInput =
            new Input<>("psi", "The sampling rate, i.e. the rate of unscheduled sequenced sampling. If you want to have rates that change over time you will also need the psiChangeTimes to be set.", (RealParameter) null);

    public Input<RealParameter> psiChangeTimesInput =
            new Input<>("psiChangeTimes", "The times at which the value of psi changes. These should be given as backwards times treating the final observation *in the tree* as the present (time zero). If psi is constant then this parameter can safely be left as the default null value.", (RealParameter) null);

    public Input<RealParameter> rhoInput =
            new Input<>("rho", "The probability of sampling lineages in a scheduled sample.", (RealParameter) null);

    public Input<RealParameter> omegaInput =
            new Input<>("omega", "The occurrence rate, i.e. the rate of unscheduled unsequenced sampling. Default value of zero (no occurrence sampling). If you want to have rates that change over time you will also need the psiChangeTimes to be set.", (RealParameter) null);

    public Input<RealParameter> omegaChangeTimesInput =
            new Input<>("omegaChangeTimes", "The times at which the value of omega changes. These should be given as backwards times treating the final observation *in the tree* as the present (time zero). If omega is constant then this parameter can safely be left as the default null value.", (RealParameter) null);

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
    protected Double[] lambdaChangeTimes;
    protected Double[] muChangeTimes;
    protected Double[] psiChangeTimes;
    protected Double rho;
    protected Double[] omegaChangeTimes;
    protected Double nu;
    protected Double originTime;

    private Double[] rateChangeTimes;

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

    // START SNEAKY
    private int intIx;
    private double[] intervalStartTimes;
    private double[] intervalEndTimes;
    private double[] lncs;
    private double[] lnls;
    private double[] lambdaValues;
    private double[] muValues;
    private double[] psiValues;
    private double[] omegaValues;
    private int[] kValues;
    // END SNEAKY

    // Unclear why it is necessary, but BEAST expects there to be a zero-argument
    // constructor and if there isn't one it freaks out.
    public TimTam() {

    }

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        this.tree = (Tree) treeInput.get();

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

        this.lambdaChangeTimes =
                (lambdaChangeTimesInput.get() != null) ? lambdaChangeTimesInput.get().getValues() : new Double[]{};
        this.muChangeTimes =
                (muChangeTimesInput.get() != null) ? muChangeTimesInput.get().getValues() : new Double[]{};
        this.psiChangeTimes =
                (psiChangeTimesInput.get() != null) ? psiChangeTimesInput.get().getValues() : new Double[]{};
        this.omegaChangeTimes =
                (omegaChangeTimesInput.get() != null) ? omegaChangeTimesInput.get().getValues() : new Double[]{};
        this.rateChangeTimes = new Double[
                this.lambdaChangeTimes.length +
                this.muChangeTimes.length +
                this.psiChangeTimes.length +
                this.omegaChangeTimes.length];
        for (int ix = 0; ix < this.rateChangeTimes.length; ix++) {
            if (ix < this.lambdaChangeTimes.length) {
                this.rateChangeTimes[ix] = this.lambdaChangeTimes[ix];
            } else if (ix < this.lambdaChangeTimes.length + this.muChangeTimes.length) {
                this.rateChangeTimes[ix] = this.muChangeTimes[ix-this.lambdaChangeTimes.length];
            } else if (ix < this.lambdaChangeTimes.length + this.muChangeTimes.length + this.psiChangeTimes.length) {
                this.rateChangeTimes[ix] = this.psiChangeTimes[ix-this.lambdaChangeTimes.length-this.muChangeTimes.length];
            } else  {
                this.rateChangeTimes[ix] = this.omegaChangeTimes[ix-this.lambdaChangeTimes.length-this.muChangeTimes.length-this.psiChangeTimes.length];
            }
        }
        Arrays.sort(this.rateChangeTimes);

        // We need to be careful in counting the number of time intervals because there may be multiple
        // tree nodes counted as a single catastrophe.
        this.numTimeIntervals =
                this.rateChangeTimes.length +
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

        updateRateAndProbParams();

        // START SNEAKY
        this.lncs = new double[this.numTimeIntervals];
        this.lnls = new double[this.numTimeIntervals];
        this.lambdaValues = new double[this.numTimeIntervals];
        this.muValues = new double[this.numTimeIntervals];
        this.psiValues = new double[this.numTimeIntervals];
        this.omegaValues = new double[this.numTimeIntervals];

        this.intervalStartTimes = new double[this.numTimeIntervals];
        this.intervalEndTimes = new double[this.numTimeIntervals];
        this.kValues = new int[this.numTimeIntervals];
        this.intervalStartTimes[0] = originTime;
        this.intervalEndTimes[0] = this.intervalTerminators[0].getBwdTime();
        this.kValues[0] = 1;
        this.lambdaValues[0] = this.lambdaInput.get().getValues()[0];
        this.muValues[0] = this.muInput.get().getValues()[0];
        this.psiValues[0] = this.psiInput.get().getValues()[0];
        if (this.omegaInput.get() != null) {
            this.omegaValues[0] = this.omegaInput.get().getValues()[0];
        } else {
            this.omegaValues[0] = 0.0;
        }

        String tt;
        for (int ix = 1; ix < this.numTimeIntervals; ix++) {
            this.intervalStartTimes[ix] = this.intervalTerminators[ix-1].getBwdTime();
            this.intervalEndTimes[ix] = this.intervalTerminators[ix].getBwdTime();

            this.lambdaValues[ix] = this.birth(this.intervalStartTimes[ix]);
            this.muValues[ix] = this.death(this.intervalStartTimes[ix]);
            this.psiValues[ix] = this.psi(this.intervalStartTimes[ix]);
            this.omegaValues[ix] = this.omega(this.intervalStartTimes[ix]);

            tt = this.intervalTerminators[ix-1].getType();
            if (Objects.equals(tt, "birth")) {
                this.kValues[ix] = this.kValues[ix-1] + 1;
            } else if (Objects.equals(tt, "sample")) {
                this.kValues[ix] = this.kValues[ix-1] - 1;
            } else if (
                    Objects.equals(tt, "occurrence") |
                    Objects.equals(tt, "rateChange") |
                    Objects.equals(tt, "disaster")
            ) {
                this.kValues[ix] = this.kValues[ix-1];
            } else if (Objects.equals(tt, "catastrophe")) {
                this.kValues[ix] = this.kValues[ix-1] - this.intervalTerminators[ix].getCount();
            } else {
                throw new RuntimeException("Unexpected interval terminator type: " + tt);
            }
        }
        // END SNEAKY

        this.nb.setZero();
        this.timeFromOriginToFinalDatum = this.originTime - this.intervalEndTimes[this.numTimeIntervals-1];
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

    public double birth(int intervalIx) {
        return this.lambdaValues[intervalIx];
    }

    public double birth(double bwdTime) {
        Double[] lambda = this.lambdaInput.get().getValues();
        if (lambda.length > 1) {
            double cBwdTime = Double.POSITIVE_INFINITY;
            double cLambda = lambda[0];
            double nBwdTime = this.lambdaChangeTimes[0];
            int ix = 0;
            while (true) {
                if (cBwdTime >= bwdTime & bwdTime > nBwdTime) {
                    return cLambda;
                } else {
                    ix+=1;
                    cBwdTime = this.lambdaChangeTimes[ix-1];
                    cLambda = lambda[ix];
                    if (ix < this.lambdaChangeTimes.length) {
                        nBwdTime = this.lambdaChangeTimes[ix];
                    } else {
                        nBwdTime = Double.NEGATIVE_INFINITY;
                    }
                }
            }
        } else {
            return lambda[0];
        }
    }

    public double death(int intervalIx) {
        return this.muValues[intervalIx];
    }

    public double death(double bwdTime) {
        Double[] mu = this.muInput.get().getValues();
        if (mu.length > 1) {
            double cBwdTime = Double.POSITIVE_INFINITY;
            double cMu = mu[0];
            double nBwdTime = this.muChangeTimes[0];
            int ix = 0;
            while (true) {
                if (cBwdTime >= bwdTime & bwdTime > nBwdTime) {
                    return cMu;
                } else {
                    ix+=1;
                    cBwdTime = this.muChangeTimes[ix-1];
                    cMu = mu[ix];
                    if (ix < this.muChangeTimes.length) {
                        nBwdTime = this.muChangeTimes[ix];
                    } else {
                        nBwdTime = Double.NEGATIVE_INFINITY;
                    }
                }
            }
        } else {
            return mu[0];
        }
    }

    public double psi(int intervalIx) {
        return this.psiValues[intervalIx];
    }

    public double psi(double bwdTime) {
        Double[] psi = this.psiInput.get().getValues();
        if (psi.length > 1) {
            double cBwdTime = Double.POSITIVE_INFINITY;
            double cPsi = psi[0];
            double nBwdTime = this.psiChangeTimes[0];
            int ix = 0;
            while (true) {
                if (cBwdTime >= bwdTime & bwdTime > nBwdTime) {
                    return cPsi;
                } else {
                    ix+=1;
                    cBwdTime = this.psiChangeTimes[ix-1];
                    cPsi = psi[ix];
                    if (ix < this.psiChangeTimes.length) {
                        nBwdTime = this.psiChangeTimes[ix];
                    } else {
                        nBwdTime = Double.NEGATIVE_INFINITY;
                    }
                }
            }
        } else {
            return psi[0];
        }
    }

    public double omega(int intervalIx) {
        return this.omegaValues[intervalIx];
    }

    public double omega(double bwdTime) {
        if (this.omegaInput.get() == null) {
            return 0.0;
        }
        Double[] omega = this.omegaInput.get().getValues();
        if (omega.length > 1) {
            double cBwdTime = Double.POSITIVE_INFINITY;
            double cOmega = omega[0];
            double nBwdTime = this.omegaChangeTimes[0];
            int ix = 0;
            while (true) {
                if (cBwdTime >= bwdTime & bwdTime > nBwdTime) {
                    return cOmega;
                } else {
                    ix+=1;
                    cBwdTime = this.omegaChangeTimes[ix-1];
                    cOmega = omega[ix];
                    if (ix < this.omegaChangeTimes.length) {
                        nBwdTime = this.omegaChangeTimes[ix];
                    } else {
                        nBwdTime = Double.NEGATIVE_INFINITY;
                    }
                }
            }
        } else {
            return omega[0];
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
        updateRateAndProbParams();
        return calculateTreeLogLikelihood();
    }

    /**
     * Update the (primitive) local parameters to match the input.
     */
    private void updateRateAndProbParams() {
        this.rho = (rhoInput.get() != null) ? rhoInput.get().getValue() : null;
        this.nu = (nuInput.get() != null) ? nuInput.get().getValue() : null;

    }

    // SNEAKY
    public double sneakySum(double[] xs) {
        double tmp = 0;
        for (int ix = 0; ix < xs.length; ix++) {
            tmp += xs[ix];
        }
        return(tmp);
    }

    /**
     * Generic likelihood calculation
     */
    public final double calculateTreeLogLikelihood() {
        this.nb.setZero();

        this.intIx = 0; // SNEAKY
        for (int ix = 0; ix < this.numTimeIntervals; ix++) {
            processInterval(
                    this.intIx,
                    this.kValues[ix]
            );
            processObservation(this.intervalTerminators[ix], this.kValues[ix]);
            this.intIx++; // SNEAKY
        }

        // if the likelihood conditions upon the observation of the process then we need to account for this in the
        // log-likelihood.
        if (this.conditionOnObservation) {
            double probUnobserved = p0(this.timeFromOriginToFinalDatum, 0.0, 1.0);
            return sneakySum(this.lnls) + sneakySum(this.lncs) - Math.log(1 - probUnobserved);
        } else {
            return sneakySum(this.lnls) + sneakySum(this.lncs);
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
     * @param k the number of lineages in the reconstructed tree just prior to the observation.
     *
     * @see TimTamIntervalTerminator
     *
     */
    private void processObservation(TimTamIntervalTerminator intTerminator, int k) {
        String intTypeStr = intTerminator.getType();
        Double obsBwdTime = intTerminator.getBwdTime();

        if (Objects.equals(intTypeStr, "birth")) {
            this.lnls[this.intIx] = Math.log(birth(obsBwdTime));
        } else if (Objects.equals(intTypeStr, "sample")) {
            this.lnls[this.intIx] = Math.log(psi(obsBwdTime));
        } else if (Objects.equals(intTypeStr, "occurrence")) {
            this.lnls[this.intIx] = Math.log(omega(obsBwdTime)) + this.nb.lnPGFDash1(1.0);
            double lnRp1 = Math.log(Math.exp(this.nb.getLnR()) + 1.0);
            this.nb.setLnPAndLnR(this.nb.getLnP(), lnRp1);
        } else if (Objects.equals(intTypeStr, "catastrophe")) {
            int n = intTerminator.getCount();
            double rho = rho();
            this.lnls[this.intIx] = (k - n) * Math.log(1 - rho)
                    + n * Math.log(rho)
                    + this.nb.lnPGF(1 - rho);
            this.nb.setLnPAndLnR(Math.log(1 - rho) + this.nb.getLnP(),
                    this.nb.getLnR());
        } else if (Objects.equals(intTypeStr, "disaster")) {
            int h = intTerminator.getCount();
            double nu = nu();
            if (h > 0) {
                this.lnls[this.intIx] = k * Math.log(1 - nu)
                        + h * Math.log(nu)
                        + this.nb.lnPGFDash(h, 1 - nu);
            } else if (h == 0) {
                this.lnls[this.intIx] = k * Math.log(1 - nu)
                        + this.nb.lnPGFDash(h, 1 - nu); // this should just be the lnPGF.
            } else {
                throw new RuntimeException("a disaster with a negative number of cases should never happen.");
            }
            this.nb.setLnPAndLnR(
                    Math.log(1 - nu) + this.nb.getLnP(),
                    Math.log(Math.exp(this.nb.getLnR()) + h));
        } else if (Objects.equals(intTypeStr, "rateChange")) {
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
     * @param bwdTimeIntervalEnd the time at which the interval ended
     * @param k the number of lineages in the reconstructed tree during the interval.
     */
//    private void processInterval(double intervalDuration, double bwdTimeIntervalEnd, int k) {
    private void processInterval(int intervalIx, int k) {

        double lnC, lnMean, lnVariance;

        double p0Val = p0(intervalIx, 1.0);
        double lnP0Dash1Val = lnP0Dash1(intervalIx);
        double lnP0Dash2Val = lnP0Dash2(intervalIx);
        double lnRVal = lnR(intervalIx);
        double lnRDash1Val = lnRDash1(intervalIx);
        double lnRDash2Val = lnRDash2(intervalIx);

        double lnFM0, lnFM1, lnFM2;
        if (!this.nb.isZero) {
            assert k >= 0;
            if (k > 0) {

                lnFM0 = this.nb.lnPGF(p0Val) + k * lnRVal;
                tmpArry[0] = this.nb.lnPGFDash1(p0Val) + lnP0Dash1Val + k * lnRVal;
                tmpArry[1] = Math.log(k) + (k-1) * lnRVal + lnRDash1Val + this.nb.lnPGF(p0Val);
                lnFM1 = logSumExp(tmpArry, 2);
                tmpArry[0] = this.nb.lnPGFDash2(p0Val) + 2 * lnP0Dash1Val + k * lnRVal;
                tmpArry[1] = this.nb.lnPGFDash1(p0Val) + lnP0Dash2Val + k * lnRVal;
                tmpArry[2] = Math.log(2) + this.nb.lnPGFDash1(p0Val) + lnP0Dash1Val + Math.log(k) + (k - 1) * lnRVal + lnRDash1Val;
                tmpArry[3] = this.nb.lnPGF(p0Val) + Math.log(k) + Math.log(k - 1) + (k - 2) * lnRVal + 2 * lnRDash1Val;
                tmpArry[4] = this.nb.lnPGF(p0Val) + Math.log(k) + (k-1) * lnRVal + lnRDash2Val;
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

        this.lncs[this.intIx] = lnC; // SNEAKY
    }

    private double lnR(int intervalIx) {
        odeHelpers(intervalIx);
        return Math.log(this.ohDiscriminant)
                + Math.log(this.ohExpFact)
                - (2 * Math.log(birth(intervalIx)))
                - (2 * Math.log((this.ohX2 - 1.0) - (this.ohX1 - 1.0) * this.ohExpFact));
    }

    private double lnRDash1(int intervalIx) {
        odeHelpers(intervalIx);
        return Math.log(2)
                + Math.log(1 - this.ohExpFact)
                + Math.log(this.ohExpFact)
                + Math.log(this.ohDiscriminant)
                - 2 * Math.log(birth(intervalIx))
                - 3 * Math.log((this.ohX2 - this.ohExpFact * (this.ohX1 - 1.0) - 1.0));
    }

    private double lnRDash2(int intervalIx) {
        odeHelpers(intervalIx);
        return Math.log(6)
                + 2 * Math.log(1 - this.ohExpFact)
                + Math.log(this.ohExpFact)
                + Math.log(this.ohDiscriminant)
                - 2 * Math.log(birth(intervalIx))
                - 4 * Math.log((this.ohX2 - this.ohExpFact * (this.ohX1 - 1.0) - 1.0));
    }

    protected double p0(double intervalDuration, double bwdTimeIntervalEnd, double z) {
        odeHelpers(intervalDuration, bwdTimeIntervalEnd);
        return (this.ohX1 * (this.ohX2 - z) - this.ohX2 * (this.ohX1 - z) * this.ohExpFact) / ((this.ohX2 - z) - (this.ohX1 - z) * this.ohExpFact);
    }

    protected double p0(int intervalIx, double z) {
        odeHelpers(intervalIx);
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
    protected double lnP0Dash1(int intervalIx) {
        odeHelpers(intervalIx);
        double aa = this.ohX2 - this.ohX1 * this.ohExpFact;
        double bb = 1 - this.ohExpFact;
        double cc = this.ohX2 * this.ohExpFact - this.ohX1;
        return Math.log(cc * aa + this.ohX1 * this.ohX2 * Math.pow(bb, 2.0))
                - 2.0 * Math.log(aa - bb);
    }

    private double lnP0Dash2(int intervalIx) {
        odeHelpers(intervalIx);
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
        double bwdRateTime = bwdTimeIntervalEnd + 0.5 * intervalDuration;
        // the 0.5*intervalDuration here is to ensure that this takes the value in the middle of the interval.
        double gamma = birth(bwdRateTime) + death(bwdRateTime) + psi(bwdRateTime) + omega(bwdRateTime);
        this.ohDiscriminant = Math.pow(gamma, 2.0) - 4.0 * birth(bwdRateTime) * death(bwdRateTime);
        double sqrtDisc = Math.sqrt(this.ohDiscriminant);
        this.ohX1 = (gamma - sqrtDisc) / (2 * birth(bwdRateTime));
        this.ohX2 = (gamma + sqrtDisc) / (2 * birth(bwdRateTime));
        this.ohExpFact = Math.exp(- sqrtDisc * intervalDuration);
    }

    private void odeHelpers(int intervalIx) {
        double intervalDuration = this.intervalStartTimes[intervalIx] - this.intervalEndTimes[intervalIx];
        double gamma = birth(intervalIx) + death(intervalIx) + psi(intervalIx) + omega(intervalIx);
        this.ohDiscriminant = Math.pow(gamma, 2.0) - 4.0 * birth(intervalIx) * death(intervalIx);
        double sqrtDisc = Math.sqrt(this.ohDiscriminant);
        this.ohX1 = (gamma - sqrtDisc) / (2 * birth(intervalIx));
        this.ohX2 = (gamma + sqrtDisc) / (2 * birth(intervalIx));
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
