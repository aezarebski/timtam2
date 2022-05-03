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
 * <p>Tree prior for birth-death-sampling while tracking the distribution of hidden
 * lineages. This used BirthDeathSerialSampling.java as a starting point. There
 * is a nested class, NegativeBinomial, which is the approximation used for the
 * number of unobserved lineages in the likelihood.</p>
 *
 * <p>Time is measured backwards from the last sequenced sample (which is treated
 * as having been collected at time zero). This is important because the times of
 * any unsequenced observations or changes in parameters need their times provided
 * this way.</p>
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
            new Input<>("rho", "The probability of sampling lineages in a scheduled sample, i.e. the probability of an individual being removed and sequenced in a scheduled sample. If you want to have this probability change through time you will also need to the rhoChangeTimes to be set.", (RealParameter) null);

    public Input<RealParameter> rhoChangeTimesInput =
            new Input<>("rhoChangeTimes", "The times at which the value of rho changes. These should be given as backwards times treating the final observation *in the tree* as the present (time zero). If rho is constant then this parameter can safely be left as the default null value.", (RealParameter) null);

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
            new Input<>("nu", "the probability of unsequenced scheduled sampling, i.e. the probability of an individual being removed but not sequenced in a scheduled sample. If you want to have this probability change through time you will also need to the nuChangeTimes to be set.", Input.Validate.OPTIONAL);

    public Input<RealParameter> nuChangeTimesInput =
            new Input<>("nuChangeTimes", "The times at which the value of nu changes. These should be given as backwards times treating the final observation *in the tree* as the present (time zero). If nu is constant then this parameter can safely be left as the default null value.", (RealParameter) null);

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
    protected Double[] rhoChangeTimes;
    protected Double[] omegaChangeTimes;
    protected Double[] nuChangeTimes;
    protected RealParameter originTime;

    private double[] paramChangeTimes;

    protected double[] catastropheTimes;
    protected int[] catastropheSizes;
    private int totalCatastropheSizes;

    // A disaster is a scheduled sample of lineages where sampled lineages are *not* sequenced so do not appear in the
    // reconstructed tree. Typically, these data will form a time series.
    protected double[] disasterTimes;
    protected int[] disasterSizes;

    // the times at which there was an occurrence sample measured in backwards time with the final tip in the tree being
    // used as zero. An occurrence is an unscheduled and unsequenced sample. These data form a point process of events.
    protected double[] occurrenceTimes;

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

    private double[] intervalStartTimes;
    private double[] intervalEndTimes;
    private double[] lncs;
    private double[] lnls;
    private double[] lambdaValues;
    private double[] muValues;
    private double[] psiValues;
    private double[] rhoValues;
    private double[] omegaValues;
    private double[] nuValues;
    private int[] kValues;

    // Unclear why it is necessary, but BEAST expects there to be a zero-argument
    // constructor and if there isn't one it freaks out.
    public TimTam() {

    }

    /**
     * <p>This function gets called once at the start of an MCMC run so any heavy
     * pre-calculations should be done here rather than in the likelihood
     * calculation.</p>
     */
    @Override
    public void initAndValidate() {
        super.initAndValidate();

        this.tree = (Tree) treeInput.get();
        this.originTime = originTimeInput.get();
        if (this.tree.getRoot().getHeight() >= this.originTime.getValue()) {
            throw new RuntimeException("tree has a root which comes before the originTime.");
        }

        if (catastropheTimesInput.get() != null) {
            this.catastropheTimes = catastropheTimesInput.get().getDoubleValues();
            this.catastropheSizes = new int[this.catastropheTimes.length];
            Node[] nodes = this.tree.getNodesAsArray();
            for (int ix = 0; ix < this.catastropheTimes.length; ix++) {
                this.catastropheSizes[ix] = 0;
                for (int jx = 0; jx < this.tree.getLeafNodeCount(); jx++) {
                    if (Math.abs(nodes[jx].getHeight() - this.catastropheTimes[ix]) < this.timeEpsilon) {
                        this.catastropheSizes[ix]++;
                    }
                }
            }
            this.totalCatastropheSizes = arraySum(this.catastropheSizes);
        } else {
            this.catastropheTimes = new double[] {};
            this.catastropheSizes = new int[] {};
            this.totalCatastropheSizes = 0;
        }

        this.occurrenceTimes = occurrenceTimesInput.get() != null ? occurrenceTimesInput.get().getDoubleValues() : new double[] {};

        if (disasterTimesInput.get() != null) {
            this.disasterTimes = disasterTimesInput.get().getDoubleValues();
            this.disasterSizes = new int[disasterSizesInput.get().getValues().length];
            for (int ix = 0; ix < this.disasterSizes.length; ix++) {
                this.disasterSizes[ix] = disasterSizesInput.get().getNativeValue(ix);
            }
        } else {
            this.disasterTimes = new double[] {};
            this.disasterSizes = new int[] {};
        }

        this.conditionOnObservation = conditionOnObservationInput.get();

        this.lambdaChangeTimes =
                (lambdaChangeTimesInput.get() != null) ? lambdaChangeTimesInput.get().getValues() : new Double[]{};
        this.muChangeTimes =
                (muChangeTimesInput.get() != null) ? muChangeTimesInput.get().getValues() : new Double[]{};
        this.psiChangeTimes =
                (psiChangeTimesInput.get() != null) ? psiChangeTimesInput.get().getValues() : new Double[]{};
        this.rhoChangeTimes =
                (rhoChangeTimesInput.get() != null) ? rhoChangeTimesInput.get().getValues() : new Double[]{};
        this.omegaChangeTimes =
                (omegaChangeTimesInput.get() != null) ? omegaChangeTimesInput.get().getValues() : new Double[]{};
        this.nuChangeTimes =
                (nuChangeTimesInput.get() != null) ? nuChangeTimesInput.get().getValues() : new Double[]{};
        this.paramChangeTimes = new double[
                this.lambdaChangeTimes.length +
                this.muChangeTimes.length +
                this.psiChangeTimes.length +
                this.rhoChangeTimes.length +
                this.omegaChangeTimes.length +
                this.nuChangeTimes.length];
        for (int ix = 0; ix < this.paramChangeTimes.length; ix++) {
            if (ix < this.lambdaChangeTimes.length) {
                this.paramChangeTimes[ix] = this.lambdaChangeTimes[ix];
            } else if (ix < this.lambdaChangeTimes.length + this.muChangeTimes.length) {
                this.paramChangeTimes[ix] = this.muChangeTimes[ix-this.lambdaChangeTimes.length];
            } else if (ix < this.lambdaChangeTimes.length + this.muChangeTimes.length + this.psiChangeTimes.length) {
                this.paramChangeTimes[ix] = this.psiChangeTimes[ix-this.lambdaChangeTimes.length-this.muChangeTimes.length];
            } else if (ix < this.lambdaChangeTimes.length + this.muChangeTimes.length + this.psiChangeTimes.length + this.rhoChangeTimes.length) {
                this.paramChangeTimes[ix] = this.rhoChangeTimes[ix-this.lambdaChangeTimes.length-this.muChangeTimes.length-this.psiChangeTimes.length];
            } else if (ix < this.lambdaChangeTimes.length + this.muChangeTimes.length + this.psiChangeTimes.length + this.rhoChangeTimes.length + this.omegaChangeTimes.length) {
                this.paramChangeTimes[ix] = this.omegaChangeTimes[ix-this.lambdaChangeTimes.length-this.muChangeTimes.length-this.psiChangeTimes.length-this.rhoChangeTimes.length];
            } else {
                this.paramChangeTimes[ix] = this.nuChangeTimes[ix-this.lambdaChangeTimes.length-this.muChangeTimes.length-this.psiChangeTimes.length-this.rhoChangeTimes.length-this.omegaChangeTimes.length];
            }
        }
        Arrays.sort(this.paramChangeTimes);

        // We need to be careful in counting the number of time intervals because there may be multiple
        // tree nodes counted as a single catastrophe.
        this.numTimeIntervals =
                this.paramChangeTimes.length +
                this.disasterTimes.length + this.occurrenceTimes.length +
                this.tree.getNodesAsArray().length - this.totalCatastropheSizes + this.catastropheTimes.length;
        this.intervalTerminators = new TimTamIntervalTerminator[this.numTimeIntervals];
        for (int ix = 0; ix < this.numTimeIntervals; ix++) {
            this.intervalTerminators[ix] = new TimTamIntervalTerminator();
        }
        this.intervalStartTimes = new double[this.numTimeIntervals];
        this.intervalEndTimes = new double[this.numTimeIntervals];
        updateIntervalTerminators();

        this.lncs = new double[this.numTimeIntervals];
        this.lnls = new double[this.numTimeIntervals];
        this.lambdaValues = new double[this.numTimeIntervals];
        this.muValues = new double[this.numTimeIntervals];
        this.psiValues = new double[this.numTimeIntervals];
        this.rhoValues = new double[this.numTimeIntervals];
        this.omegaValues = new double[this.numTimeIntervals];
        this.nuValues = new double[this.numTimeIntervals];
        this.kValues = new int[this.numTimeIntervals];
        updateRateAndProbParams();
    }

    /**
     * <p>Predicate for being an unscheduled node.</p>
     * <p>Because there is the possibility of numerical issues we use the timeEpsilon
     * to check if a node is sufficiently close to a catastrophe to be
     * considered one.</p>
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

    public double rho(int intervalIx) {
        return this.rhoValues[intervalIx];
    }

    public double rho(double bwdTime) {
        if (this.rhoInput.get() == null) {
            return 0.0;
        }
        Double[] rho = this.rhoInput.get().getValues();
        if (rho.length > 1) {
            double cBwdTime = Double.POSITIVE_INFINITY;
            double cRho = rho[0];
            double nBwdTime = this.rhoChangeTimes[0];
            int ix = 0;
            while (true) {
                if (cBwdTime >= bwdTime & bwdTime > nBwdTime) {
                    return cRho;
                } else {
                    ix+=1;
                    cBwdTime = this.rhoChangeTimes[ix-1];
                    cRho = rho[ix];
                    if (ix < this.rhoChangeTimes.length) {
                        nBwdTime = this.rhoChangeTimes[ix];
                    } else {
                        nBwdTime = Double.NEGATIVE_INFINITY;
                    }
                }
            }
        } else {
            return rho[0];
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

    public double nu(int intervalIx) {
        return this.nuValues[intervalIx];
    }

    public double nu(double bwdTime) {
        if (this.nuInput.get() == null) {
            return 0.0;
        }
        Double[] nu = this.nuInput.get().getValues();
        if (nu.length > 1) {
            double cBwdTime = Double.POSITIVE_INFINITY;
            double cNu = nu[0];
            double nBwdTime = this.nuChangeTimes[0];
            int ix = 0;
            while (true) {
                if (cBwdTime >= bwdTime & bwdTime > nBwdTime) {
                    return cNu;
                } else {
                    ix+=1;
                    cBwdTime = this.nuChangeTimes[ix-1];
                    cNu = nu[ix];
                    if (ix < this.nuChangeTimes.length) {
                        nBwdTime = this.nuChangeTimes[ix];
                    } else {
                        nBwdTime = Double.NEGATIVE_INFINITY;
                    }
                }
            }
        } else {
            return nu[0];
        }
    }

    @Override
    public double calculateLogP() {
        // If the tree is tall enough that the root happens before the origin then this point in parameter space has
        // probability zero.
        if (this.tree.getRoot().getHeight() >= this.originTime.getValue()) {
            return Double.NEGATIVE_INFINITY;
        }

        this.nb.setIsZero(true);
        updateIntervalTerminators();
        updateRateAndProbParams();
        for (int ix = 0; ix < this.numTimeIntervals; ix++) {
            updateOdeHelpers(ix);
            processInterval(ix);
            processObservation(ix);
        }

        // if the likelihood conditions upon the observation of the process then we need to account for this in the
        // log-likelihood.
        if (this.conditionOnObservation) {
            if (this.paramChangeTimes.length == 0) {
            double probUnobserved = p0(this.timeFromOriginToFinalDatum, 0.0, 1.0);
            return arraySum(this.lnls) + arraySum(this.lncs) - Math.log(1 - probUnobserved);
            } else {
                throw new RuntimeException("conditioning upon observation is not yet implemented for varying rates.");
            }
        } else {
            return arraySum(this.lnls) + arraySum(this.lncs);
        }
    }

    /**
     * This function updates the interval terminators which could change if there are changes to the Tree object.
     *
     * <p>The use of {@link java.util.Arrays#sort} here is appropriate according to the
     * <a href="https://docs.oracle.com/javase/7/docs/api/java/util/Arrays.html#sort(java.lang.Object[])">documentation</a>
     * because it is suitable for sorting concatenated sorted arrays.</p>
     */
    private void updateIntervalTerminators() {
        int iTx = 0;

        for (Double paramChangeTime : this.paramChangeTimes) {
            this.intervalTerminators[iTx].setTypeTimeAndCount("paramValueChange", paramChangeTime, OptionalInt.empty());
            iTx++;
        }

        for (Double occurrenceTime : this.occurrenceTimes) {
            this.intervalTerminators[iTx].setTypeTimeAndCount("occurrence", occurrenceTime, OptionalInt.empty());
            iTx++;
        }

        for (int ix = 0; ix < this.catastropheTimes.length; ix++) {
            this.intervalTerminators[iTx].setTypeTimeAndCount("catastrophe", this.catastropheTimes[ix], OptionalInt.of(this.catastropheSizes[ix]));
            iTx++;
        }

        for (int ix = 0; ix < this.disasterTimes.length; ix++) {
            this.intervalTerminators[iTx].setTypeTimeAndCount("disaster", this.disasterTimes[ix], OptionalInt.of(this.disasterSizes[ix]));
            iTx++;
        }

        for (Node node : this.tree.getNodesAsArray()) {
            if (isUnscheduledTreeNode(node)) {
                this.intervalTerminators[iTx].setTypeTimeAndCount(node.isLeaf() ? "sample" : "birth", node.getHeight(), OptionalInt.empty());
                iTx++;
            }
        }
        Arrays.sort(this.intervalTerminators);
        this.intervalStartTimes[0] = this.originTime.getValue();
        this.intervalEndTimes[0] = this.intervalTerminators[0].getBwdTime();
        for (int ix = 1; ix < this.numTimeIntervals; ix++) {
            this.intervalStartTimes[ix] = this.intervalTerminators[ix-1].getBwdTime();
            this.intervalEndTimes[ix] = this.intervalTerminators[ix].getBwdTime();
        }
        this.timeFromOriginToFinalDatum = this.originTime.getValue() - this.intervalEndTimes[this.numTimeIntervals-1];
    }

    /**
     * Update the (primitive) local parameters to match the input.
     */
    private void updateRateAndProbParams() {
        this.kValues[0] = 1;
        this.lambdaValues[0] = this.lambdaInput.get().getValues()[0];
        this.muValues[0] = this.muInput.get().getValues()[0];
        this.psiValues[0] = this.psiInput.get().getValues()[0];
        if (this.rhoInput.get() != null) {
            this.rhoValues[0] = this.rhoInput.get().getValues()[0];
        } else {
            this.rhoValues[0] = 0.0;
        }
        if (this.omegaInput.get() != null) {
            this.omegaValues[0] = this.omegaInput.get().getValues()[0];
        } else {
            this.omegaValues[0] = 0.0;
        }
        if (this.nuInput.get() != null) {
            this.nuValues[0] = this.nuInput.get().getValues()[0];
        } else {
            this.nuValues[0] = 0.0;
        }

        String tt;
        for (int ix = 1; ix < this.numTimeIntervals; ix++) {
            this.lambdaValues[ix] = this.birth(this.intervalStartTimes[ix]);
            this.muValues[ix] = this.death(this.intervalStartTimes[ix]);
            this.psiValues[ix] = this.psi(this.intervalStartTimes[ix]);
            this.rhoValues[ix] = this.rho(this.intervalStartTimes[ix]);
            this.omegaValues[ix] = this.omega(this.intervalStartTimes[ix]);
            this.nuValues[ix] = this.nu(this.intervalStartTimes[ix]);

            tt = this.intervalTerminators[ix-1].getType();
            if (Objects.equals(tt, "birth")) {
                this.kValues[ix] = this.kValues[ix-1] + 1;
            } else if (Objects.equals(tt, "sample")) {
                this.kValues[ix] = this.kValues[ix-1] - 1;
            } else if (
                    Objects.equals(tt, "occurrence") |
                            Objects.equals(tt, "paramValueChange") |
                            Objects.equals(tt, "disaster")
            ) {
                this.kValues[ix] = this.kValues[ix-1];
            } else if (Objects.equals(tt, "catastrophe")) {
                this.kValues[ix] = this.kValues[ix-1] - this.intervalTerminators[ix-1].getCount();
            } else if (Objects.equals(tt, "disaster")) {
                this.kValues[ix] = this.kValues[ix-1] - this.intervalTerminators[ix-1].getCount();
            } else {
                throw new RuntimeException("Unexpected interval terminator type: " + tt);
            }
        }
    }

    public double arraySum(double[] xs) {
        double tmp = 0.0;
        for (int ix = 0; ix < xs.length; ix++) {
            tmp += xs[ix];
        }
        return(tmp);
    }

    public int arraySum(int[] xs) {
        int tmp = 0;
        for (int ix = 0; ix < xs.length; ix++) {
            tmp += xs[ix];
        }
        return(tmp);
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
    private void processObservation(int intNum) {
        TimTamIntervalTerminator intTerminator = this.intervalTerminators[intNum];
        int k = this.kValues[intNum];
        String intTypeStr = intTerminator.getType();

        if (Objects.equals(intTypeStr, "birth")) {
            this.lnls[intNum] = Math.log(birth(intNum));
        } else if (Objects.equals(intTypeStr, "sample")) {
            this.lnls[intNum] = Math.log(psi(intNum));
        } else if (Objects.equals(intTypeStr, "occurrence")) {
            this.lnls[intNum] = Math.log(omega(intNum)) + this.nb.lnPGFDash1(1.0);
            double lnRp1 = Math.log(Math.exp(this.nb.getLnR()) + 1.0);
            this.nb.setLnPAndLnR(this.nb.getLnP(), lnRp1);
        } else if (Objects.equals(intTypeStr, "catastrophe")) {
            int n = intTerminator.getCount();
            this.lnls[intNum] = (k - n) * Math.log(1 - rho(intNum))
                    + n * Math.log(rho(intNum))
                    + this.nb.lnPGF(1 - rho(intNum));
            this.nb.setLnPAndLnR(Math.log(1 - rho(intNum)) + this.nb.getLnP(),
                    this.nb.getLnR());
        } else if (Objects.equals(intTypeStr, "disaster")) {
            int h = intTerminator.getCount();
            if (h > 0) {
                this.lnls[intNum] = k * Math.log(1 - nu(intNum))
                        + h * Math.log(nu(intNum))
                        + this.nb.lnPGFDash(h, 1 - nu(intNum));
            } else if (h == 0) {
                this.lnls[intNum] = k * Math.log(1 - nu(intNum))
                        + this.nb.lnPGFDash(h, 1 - nu(intNum)); // this should just be the lnPGF.
            } else {
                throw new RuntimeException("a disaster with a negative number of cases should never happen.");
            }
            this.nb.setLnPAndLnR(
                    Math.log(1 - nu(intNum)) + this.nb.getLnP(),
                    Math.log(Math.exp(this.nb.getLnR()) + h));
        } else if (Objects.equals(intTypeStr, "paramValueChange")) {
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
    private void processInterval(int intNum) {
        int k = this.kValues[intNum];
        double lnC, lnMean, lnVariance;
        double p0Val = p0(intNum, 1.0);
        double lnP0Dash1Val = lnP0Dash1(intNum);
        double lnP0Dash2Val = lnP0Dash2(intNum);
        double lnRVal = lnR(intNum);
        double lnRDash1Val = lnRDash1(intNum);
        double lnRDash2Val = lnRDash2(intNum);

        double lnPGFVal = this.nb.lnPGF(p0Val);
        double lnPGFDash1Val = this.nb.lnPGFDash1(p0Val);
        double lnPGFDash2Val = this.nb.lnPGFDash2(p0Val);

        double lnFM0, lnFM1, lnFM2;
        if (!this.nb.getIsZero()) {
            assert k >= 0;
            if (k > 0) {
                lnFM0 = lnPGFVal + k * lnRVal;
                tmpArry[0] = lnPGFDash1Val + lnP0Dash1Val + k * lnRVal;
                tmpArry[1] = Math.log(k) + (k-1) * lnRVal + lnRDash1Val + lnPGFVal;
                lnFM1 = logSumExp(tmpArry, 2);
                tmpArry[0] = lnPGFDash2Val + 2 * lnP0Dash1Val + k * lnRVal;
                tmpArry[1] = lnPGFDash1Val + lnP0Dash2Val + k * lnRVal;
                tmpArry[2] = Math.log(2) + lnPGFDash1Val + lnP0Dash1Val + Math.log(k) + (k - 1) * lnRVal + lnRDash1Val;
                tmpArry[3] = lnPGFVal + Math.log(k) + Math.log(k - 1) + (k - 2) * lnRVal + 2 * lnRDash1Val;
                tmpArry[4] = lnPGFVal + Math.log(k) + (k-1) * lnRVal + lnRDash2Val;
                lnFM2 = logSumExp(tmpArry, 5);
            } else {
                lnFM0 = lnPGFVal;
                lnFM1 = lnPGFDash1Val + lnP0Dash1Val;
                tmpArry[0] = lnPGFDash2Val + 2 * lnP0Dash1Val;
                tmpArry[1] = lnPGFDash1Val + lnP0Dash2Val;
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
        this.lncs[intNum] = lnC;
    }

    private double lnR(int intervalIx) {
        return Math.log(this.ohDiscriminant)
                + Math.log(this.ohExpFact)
                - (2 * Math.log(birth(intervalIx)))
                - (2 * Math.log((this.ohX2 - 1.0) - (this.ohX1 - 1.0) * this.ohExpFact));
    }

    private double lnRDash1(int intervalIx) {
        return Math.log(2)
                + Math.log(1 - this.ohExpFact)
                + Math.log(this.ohExpFact)
                + Math.log(this.ohDiscriminant)
                - 2 * Math.log(birth(intervalIx))
                - 3 * Math.log((this.ohX2 - this.ohExpFact * (this.ohX1 - 1.0) - 1.0));
    }

    private double lnRDash2(int intervalIx) {
        return Math.log(6)
                + 2 * Math.log(1 - this.ohExpFact)
                + Math.log(this.ohExpFact)
                + Math.log(this.ohDiscriminant)
                - 2 * Math.log(birth(intervalIx))
                - 4 * Math.log((this.ohX2 - this.ohExpFact * (this.ohX1 - 1.0) - 1.0));
    }

    protected double p0(double intervalDuration, double bwdTimeIntervalEnd, double z) {
        updateOdeHelpers(intervalDuration, bwdTimeIntervalEnd);
        return (this.ohX1 * (this.ohX2 - z) - this.ohX2 * (this.ohX1 - z) * this.ohExpFact) / ((this.ohX2 - z) - (this.ohX1 - z) * this.ohExpFact);
    }

    protected double p0(int intervalIx, double z) {
        return (this.ohX1 * (this.ohX2 - z) - this.ohX2 * (this.ohX1 - z) * this.ohExpFact) / ((this.ohX2 - z) - (this.ohX1 - z) * this.ohExpFact);
    }
    protected double lnP0Dash1(double intervalDuration, double bwdTimeIntervalEnd) {
        updateOdeHelpers(intervalDuration, bwdTimeIntervalEnd);
        double aa = this.ohX2 - this.ohX1 * this.ohExpFact;
        double bb = 1 - this.ohExpFact;
        double cc = this.ohX2 * this.ohExpFact - this.ohX1;
        return Math.log(cc * aa + this.ohX1 * this.ohX2 * Math.pow(bb, 2.0))
                - 2.0 * Math.log(aa - bb);
    }
    protected double lnP0Dash1(int intervalIx) {
        double aa = this.ohX2 - this.ohX1 * this.ohExpFact;
        double bb = 1 - this.ohExpFact;
        double cc = this.ohX2 * this.ohExpFact - this.ohX1;
        return Math.log(cc * aa + this.ohX1 * this.ohX2 * Math.pow(bb, 2.0))
                - 2.0 * Math.log(aa - bb);
    }

    private double lnP0Dash2(int intervalIx) {
        double aa = this.ohX2 - this.ohX1 * this.ohExpFact;
        double bb = 1 - this.ohExpFact;
        double cc = this.ohX2 * this.ohExpFact - this.ohX1;
        return Math.log(2.0)
                + Math.log(bb)
                + Math.log(cc * aa + this.ohX1 * this.ohX2 * Math.pow(bb, 2.0))
                - 3.0 * Math.log(aa - bb);
    }

    private double ohX1, ohX2, ohDiscriminant, ohExpFact;

    private void updateOdeHelpers(double intervalDuration, double bwdTimeIntervalEnd) {
        double bwdRateTime = bwdTimeIntervalEnd + 0.5 * intervalDuration;
        // the 0.5*intervalDuration here is to ensure that this takes the value in the middle of the interval.
        double gamma = birth(bwdRateTime) + death(bwdRateTime) + psi(bwdRateTime) + omega(bwdRateTime);
        this.ohDiscriminant = Math.pow(gamma, 2.0) - 4.0 * birth(bwdRateTime) * death(bwdRateTime);
        double sqrtDisc = Math.sqrt(this.ohDiscriminant);
        this.ohX1 = (gamma - sqrtDisc) / (2 * birth(bwdRateTime));
        this.ohX2 = (gamma + sqrtDisc) / (2 * birth(bwdRateTime));
        this.ohExpFact = Math.exp(- sqrtDisc * intervalDuration);
    }

    /**
     * Update some expressions that are useful in multiple calculations. This function gets called once before each
     * interval is processed since the results do not change within an interval.
     *
     * @implNote The local variables <code>br</code> and <code>dr</code> are there to avoid repeated calls to the
     * function that looks these values up. The use of {@link java.lang.Math#pow} seems to be quicker than using
     * <code>gamma * gamma</code>, possibly because this is using native code.
     *
     * @param intervalIx
     */
    private void updateOdeHelpers(int intervalIx) {
        double intervalDuration = this.intervalStartTimes[intervalIx] - this.intervalEndTimes[intervalIx];
        double br = this.lambdaValues[intervalIx];
        double dr = this.muValues[intervalIx];
        double gamma = br + dr + psi(intervalIx) + omega(intervalIx);
        this.ohDiscriminant = Math.pow(gamma, 2.0) - 4.0 * br * dr;
        double sqrtDisc = Math.sqrt(this.ohDiscriminant);
        this.ohX1 = (gamma - sqrtDisc) / (2 * br);
        this.ohX2 = (gamma + sqrtDisc) / (2 * br);
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

    /**
     * This is used to report the first interval duration as a way to assess if the chain is mixing properly.
     *
     * @return the duration of the first interval.
     */
    protected double getFirstIntervalDuration() {
        return this.intervalStartTimes[0] - this.intervalEndTimes[0];
    }
}
