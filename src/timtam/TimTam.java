package timtam;


import beast.base.core.Citation;
import beast.base.core.Input;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeDistribution;

import java.util.*;

/**
 * <p>Tree prior based on the birth-death-sampling process which can also keep
 * track of the number of lineages through time. This class was originally based
 * on BirthDeathSerialSampling.java from the BDSKY package. There is a class,
 * {@link HiddenLineageDist}, which is the approximation used for the number of
 * unobserved lineages in the likelihood.</p>
 *
 * <p>There are a couple of parameterisations --- see
 * {@link TimTam#parameterisationInput} --- for this model available:</p>
 *
 * <ul>
 *     <li><code>canonical</code> (the default) is in terms of rates and
 *     probabilities</li>
 *     <li><code>r0</code> is in terms of R-naught and proportions but only
 *     applies unscheduled samples</li>
 *     <li><code>timeSeriesR0</code> is in terms of R-naught but only
 *     applies to unschedules sequences and periodic scheduled unsequenced
 *     samples (ie. a time series of cases).</li>
 * </ul>
 *
 * <p>Time is measured backwards from the last sequenced sample (which is treated
 * as having been collected at time zero). This is important because the times of
 * any unsequenced observations or changes in parameters need their times provided
 * this way.</p>
 *
 * @author Alexander E. Zarebski
 */
@Citation(value = "Zarebski AE, du Plessis L, Parag KV, Pybus OG (2022) A computationally tractable birth-death model that combines phylogenetic and epidemiological data. PLOS Computational Biology 18(2): e1009805. https://doi.org/10.1371/journal.pcbi.1009805",
        year = 2022, firstAuthorSurname = "Zarebski", DOI="10.1371/journal.pcbi.1009805")
public class TimTam extends TreeDistribution {

    public Input<RealParameter> lambdaInput =
            new Input<>("lambda",
                "The birth rate, i.e. the rate at which an infected individual spawns another " +
                    "infected individual. If you want to have rates that change over time you will " +
                    "also need the lambdaChangeTimes to be set.",
                (RealParameter) null);
    public Input<RealParameter> lambdaChangeTimesInput =
            new Input<>("lambdaChangeTimes",
                "The times at which the value of lambda changes. These should be given as " +
                    "backwards times treating the final observation *in the tree* as the present " +
                    "(time zero). If lambda is constant then this parameter can safely be left as the " +
                    "default null value.",
                (RealParameter) null);
    public Input<RealParameter> muInput =
            new Input<>("mu",
                "The death rate, i.e. the rate at which individuals are removed without being " +
                    "observed. If you want to have rates that change over time you will also need the " +
                    "muChangeTimes to be set.",
                (RealParameter) null);
    public Input<RealParameter> muChangeTimesInput =
            new Input<>("muChangeTimes",
                "The times at which the value of mu changes. These should be given as backwards " +
                    "times treating the final observation *in the tree* as the present (time zero). " +
                    "If mu is constant then this parameter can safely be left as the default null " +
                    "value.",
                (RealParameter) null);
    public Input<RealParameter> psiInput =
            new Input<>("psi",
                "The sampling rate, i.e. the rate of unscheduled sequenced sampling. If you want " +
                    "to have rates that change over time you will also need the psiChangeTimes to be " +
                    "set.",
                (RealParameter) null);
    public Input<RealParameter> psiChangeTimesInput =
            new Input<>("psiChangeTimes",
                "The times at which the value of psi changes. These should be given as backwards " +
                    "times treating the final observation *in the tree* as the present (time zero). " +
                    "If psi is constant then this parameter can safely be left as the default null " +
                    "value.",
                (RealParameter) null);
    public Input<RealParameter> omegaInput =
        new Input<>("omega",
            "The occurrence rate, i.e. the rate of unscheduled unsequenced sampling. Default " +
                "value of zero (no occurrence sampling). If you want to have rates that change " +
                "over time you will also need the psiChangeTimes to be set.",
            (RealParameter) null);
    public Input<RealParameter> omegaChangeTimesInput =
        new Input<>("omegaChangeTimes",
            "The times at which the value of omega changes. These should be given as " +
                "backwards times treating the final observation *in the tree* as the present " +
                "(time zero). If omega is constant then this parameter can safely be left as the " +
                "default null value.",
            (RealParameter) null);

    // TODO Write sensible text tips to explain these parameters.
    public Input<RealParameter> r0Input = new Input<>("r0", "PUT TEXT HERE!", (RealParameter) null);
    public Input<RealParameter> r0ChangeTimesInput = new Input<>("r0ChangeTimes", "PUT TEXT HERE!", (RealParameter) null);
    public Input<RealParameter> sigmaInput = new Input<>("sigma", "PUT TEXT HERE!", (RealParameter) null);
    public Input<RealParameter> sigmaChangeTimesInput = new Input<>("sigmaChangeTimes", "PUT TEXT HERE!", (RealParameter) null);
    public Input<RealParameter> propPsiInput = new Input<>("propPsi", "PUT TEXT HERE!", (RealParameter) null);
    public Input<RealParameter> propPsiChangeTimesInput = new Input<>("propPsiChangeTimes", "PUT TEXT HERE!", (RealParameter) null);
    public Input<RealParameter> propOmegaInput = new Input<>("propOmega", "PUT TEXT HERE!", (RealParameter) null);
    public Input<RealParameter> propOmegaChangeTimesInput = new Input<>("propOmegaChangeTimes", "PUT TEXT HERE!", (RealParameter) null);
    // These are only used for the timeSeriesR0 parameterisation
    public Input<RealParameter> propTimeSeriesInput = new Input<>("propTimeSeries", "PUT TEXT HERE!", (RealParameter) null);
    public Input<RealParameter> propTimeSeriesChangeTimesInput = new Input<>("propTimeSeriesChangeTimes", "PUT TEXT HERE!", (RealParameter) null);

    public Input<RealParameter> rhoInput =
            new Input<>("rho", "The probability of sampling lineages in a scheduled sample, i.e. the probability of an individual being removed and sequenced in a scheduled sample. If you want to have this probability change through time you will also need to the rhoChangeTimes to be set.", (RealParameter) null);
    public Input<RealParameter> rhoChangeTimesInput =
            new Input<>("rhoChangeTimes", "The times at which the value of rho changes. These should be given as backwards times treating the final observation *in the tree* as the present (time zero). If rho is constant then this parameter can safely be left as the default null value.", (RealParameter) null);
    public Input<RealParameter> nuInput =
        new Input<>("nu", "the probability of unsequenced scheduled sampling, i.e. the probability of an individual being removed but not sequenced in a scheduled sample. If you want to have this probability change through time you will also need to the nuChangeTimes to be set.", Input.Validate.OPTIONAL);
    public Input<RealParameter> nuChangeTimesInput =
        new Input<>("nuChangeTimes", "The times at which the value of nu changes. These should be given as backwards times treating the final observation *in the tree* as the present (time zero). If nu is constant then this parameter can safely be left as the default null value.", (RealParameter) null);

    public Input<RealParameter> originTimeInput =
            new Input<>("originTime",
                    "The (backwards) time of the origin (relative to the most recent observation in the tree which is considered as having occurred at time zero).",
                    (RealParameter) null);

    public Input<RealParameter> catastropheTimesInput = new Input<>("catastropheTimes", "the times at which a scheduled sequenced sample was attempted", Input.Validate.OPTIONAL);

    public Input<RealParameter> occurrenceTimesInput =
            new Input<>("occurrenceTimes",
                    "The times at which there was an occurrence (omega-sample) event. This should be entered as backwards time relative to the final time in the tree which is treated as the present (time zero). The default variable for this is null meaning no occurrence observations, leaving this empty does not mean that the omega rate is assumed to be zero.");

    public Input<RealParameter> disasterTimesInput =
            new Input<>("disasterTimes",
                    "The times at which a scheduled unsequenced sample was attempted. This should be entered as backwards time relative to the final time in the tree which is treated as the present (time zero). The default variable for this is null indicating that no scheduled unsequenced samples were attempted.",
                    (RealParameter) null);
    public Input<IntegerParameter> disasterSizesInput =
            new Input<>("disasterSizes",
                    "The size of each scheduled unsequenced sample at the corresponding disaster time.",
                    (IntegerParameter) null);

    public Input<RealParameter> historyTimesInput =
	new Input<>("historyTimes",
                    "The times at which the historical number of hidden lineages is to be estimated for." +
		    "This should be entered as backwards time relative to the final time in the tree which is treated as the present (time zero).",
                    (RealParameter) null);
    public Input<IntegerParameter> historySizesInput =
	new Input<>("historySizes",
                    "A parameter describing the number of hidden lineages at the corresponding historical time (see historyTimes parameter).",
                    (IntegerParameter) null);

    public Input<Boolean> conditionOnObservationInput =
        new Input<>("conditionOnObservation",
            "if is true then condition on sampling at least one individual (psi-sampling). " +
                "The default value is true.",
            true);

    /**
     * The name of the parameterisation of the model to be specified in the XML:
     * one of <code>canonical</code>, <code>r0</code>, or
     * <code>timeSeriesR0</code>. The default value is
     * <code>canonical</code>.
     */
    public Input<String> parameterisationInput =
        new Input<>("parameterisation",
            "the name of the parameterisation of the model, one of 'canonical', 'r0', or " +
                "'timeSeriesR0'. The default value is 'canonical'.",
            "canonical");

    /**
     * These attributes make it easier to check which parameterisation is being
     * used.
     */
    private boolean usingCanonical, usingR0, usingTimeSeriesR0;

    // we specify a threshold below which two times are considered equal.
    public final double timeEpsilon = 0.00001;

    Tree tree;

    protected String parameterisation;

    private double[] lambdaValues;
    protected Double[] lambdaChangeTimes;
    private double[] muValues;
    protected Double[] muChangeTimes;
    private double[] psiValues;
    protected Double[] psiChangeTimes;
    private double[] omegaValues;
    protected Double[] omegaChangeTimes;
    private double[] rhoValues;
    protected Double[] rhoChangeTimes;
    private double[] nuValues;
    protected Double[] nuChangeTimes;

    protected Double[] r0ChangeTimes;
    protected Double[] sigmaChangeTimes;
    protected Double[] propPsiChangeTimes;
    protected Double[] propOmegaChangeTimes;
    // This is used in the timeSeriesR0 parameterisation.
    protected Double[] propTimeSeriesChangeTimes;

    protected RealParameter originTime;

    private Double[] paramChangeTimes;

    protected double[] catastropheTimes;
    protected int[] catastropheSizes;

    // A disaster is a scheduled sample of lineages where sampled lineages are
    // *not* sequenced so do not appear in the reconstructed tree. Typically,
    // these data will form a time series.
    protected double[] disasterTimes;
    protected int[] disasterSizes;

    // when using the time series r0 parameterisation it is useful to be able to
    // refer easily to the amount of time between disasters.
    protected double disasterTimeSeriesDelta;

    // The history times and sizes are used to estimate the historical number of
    // hidden lineages.
    protected int numHistoryChecks;
    protected double[] historyTimes;
    protected int[] historySizes;

    // the times at which there was an occurrence sample measured in backwards
    // time with the final tip in the tree being used as zero. An occurrence is
    // an unscheduled and unsequenced sample. These data form a point process of
    // events.
    protected double[] occurrenceTimes;

    boolean conditionOnObservation;

    // this is used to track the distribution of hidden lineages as we process
    // the data.
    private final HiddenLineageDist nb = new HiddenLineageDist();

    // this is the time from the origin until the time of the final observation.
    private double timeFromOriginToFinalDatum;

    // this is the number of intervals of time that need to be considered when
    // evaluating the likelihood. An interval could be due to a change in a rate
    // parameter or an observation, e.g. a birth or an occurrence or a disaster.
    private int numTimeIntervals;


    // this array holds the reason each interval terminated, this could be
    // because there was an observation or a rate change.
    private IntervalTerminator[] intervalTerminators;

    private double[] intervalStartTimes;
    private double[] intervalEndTimes;
    private double[] lncs;
    private double[] lnls;
    private int[] kValues;

    /**
     * <p>This function gets called once at the start of an MCMC run so any heavy
     * pre-calculations should be done here rather than in the likelihood
     * calculation.</p>
     *
     */
    @Override
    public void initAndValidate() {
        super.initAndValidate();

        this.parameterisation = parameterisationInput.get();
        this.usingCanonical = this.parameterisation.equals("canonical");
        this.usingR0 = this.parameterisation.equals("r0");
        this.usingTimeSeriesR0 = this.parameterisation.equals("timeSeriesR0");

        if (!(this.usingCanonical | this.usingR0 | this.usingTimeSeriesR0)) {
            throw new RuntimeException(
              "The parameterisation must be one of {'canonical', 'r0', 'timeSeriesR0'}. " +
              "The value provided is '" + parameterisation + "'."
            );
        } else if (this.usingR0 &
            (this.lambdaInput.get() != null |
             this.muInput.get() != null |
             this.psiInput.get() != null |
             this.omegaInput.get() != null)) {
            throw new RuntimeException(
                "The R0 parameterisation has been used but one of the canonical parameters " +
                "({'lambda', 'mu', 'psi', 'omega'}) has a non-null value."
            );
        } else if (this.usingCanonical &
            (this.r0Input.get() != null |
             this.sigmaInput.get() != null |
             this.propPsiInput.get() != null |
             this.propOmegaInput.get() != null)) {
            throw new RuntimeException(
                "The canonical parameterisation has been used but one of the R0 parameters " +
                "({'r0', 'sigma', 'propPsi', 'propOmega'}) has a non-null value."
            );
        } else if (this.usingTimeSeriesR0 &
            (this.lambdaInput.get() != null |
             this.muInput.get() != null |
             this.psiInput.get() != null |
             this.omegaInput.get() != null |
             this.propOmegaInput.get() != null |
             this.rhoInput.get() != null |
             this.nuInput.get() != null |
             this.disasterTimesInput.get() == null |
             !disasterTimesAreUniform() |
             this.catastropheTimesInput.get() != null)) {
            throw new RuntimeException(
                "The time series R0 parameterisation has been used but something is wrong: " +
                    "one of the parameters ({'lambda', 'mu', 'psi', 'omega', 'rho', 'nu', 'propOmega'}) may have a non-null value, " +
                    "there may not be any disasters times given, " +
                    "the disaster times may not be uniformly spaced, " +
                    "there are catastrophes given."
            );
        }

        this.tree = (Tree) treeInput.get();
        this.originTime = originTimeInput.get();
        if (this.tree.getRoot().getHeight() >= this.originTime.getValue()) {
            throw new RuntimeException(
                "The tree has a root which comes before the originTime. This is an event with " +
                "probability zero so you may need to adjust your initial state to ensure this has " +
                "not happened."
            );
        }

        // We need a way to store the total number of lineages removed through
        // catastrophes because this can be tricky to get out of the tree alone.
        int totalCatastropheSizes;
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
            totalCatastropheSizes = Numerics.arraySum(this.catastropheSizes);
            if (this.catastropheSizes.length != this.catastropheTimes.length) {
                throw new RuntimeException(
                        "The number of catastrophe times given does not match the number of catastrophe sizes given."
                );
            }
        } else {
            this.catastropheTimes = new double[] {};
            this.catastropheSizes = new int[] {};
            totalCatastropheSizes = 0;
        }

        this.occurrenceTimes =
                occurrenceTimesInput.get() != null ?
                    occurrenceTimesInput.get().getDoubleValues() :
                    new double[] {};

        if (disasterTimesInput.get() != null) {
            this.disasterTimes = disasterTimesInput.get().getDoubleValues();
            // The process of storing the disaster sizes is a bit messy because
            // they are stored as unboxed values and this is not the best idea.
            this.disasterSizes = new int[disasterSizesInput.get().getValues().length];
            for (int ix = 0; ix < this.disasterSizes.length; ix++) {
                this.disasterSizes[ix] = disasterSizesInput.get().getNativeValue(ix);
            }
            if (this.disasterSizes.length != this.disasterTimes.length) {
                throw new RuntimeException(
                    "The number of disaster times given does not match the number of disaster sizes given."
                );
            }
            initDisasterTimeSeriesDelta();
        } else {
            this.disasterTimes = new double[] {};
            this.disasterSizes = new int[] {};
        }

        if (historyTimesInput.get() != null) {
            if (this.historySizesInput.get().getDimension() !=
                this.historyTimesInput.get().getDimension()) {
                throw new RuntimeException(
                    "History Sizes and History Times have different lengths."
                );
            } else {
                this.numHistoryChecks = historySizesInput.get().getValues().length;
                this.historySizes = new int[numHistoryChecks];
                this.historyTimes = new double[numHistoryChecks];
                updateHistoryChecks();
            }
        } else {
            this.numHistoryChecks = 0;
            this.historyTimes = new double[] {};
            this.historySizes = new int[] {};
        }

        this.conditionOnObservation = conditionOnObservationInput.get();

        initChangeTimes();

        // We need to be careful in counting the number of time intervals
        // because there may be multiple tree nodes counted as a single
        // catastrophe.
        this.numTimeIntervals =
                this.paramChangeTimes.length +
                this.disasterTimes.length +
                this.historyTimes.length +
                this.occurrenceTimes.length +
                this.tree.getNodesAsArray().length -
                    totalCatastropheSizes +
                    this.catastropheTimes.length;

        initIntervalTerminators();
        updateIntervalTerminators();

        initLTTArray();
        updateLTTArray();

        initRateAndProbParams();
        updateRateAndProbParams();

        this.lncs = new double[this.numTimeIntervals];
        this.lnls = new double[this.numTimeIntervals];
    }

    /**
     * Initialise the value of the times series delta after performing checks to
     * make sure this makes sense.
     */
    private void initDisasterTimeSeriesDelta() {
        if (this.usingTimeSeriesR0) {
            if (disasterTimesAreUniform()) {
                this.disasterTimeSeriesDelta =
                    this.disasterTimes[0] - this.disasterTimes[1];
            } else {
                throw new RuntimeException(
                    "Cannot use time series R0 parameterisation with the given disaster times."
                );
            }
        } else {
            this.disasterTimeSeriesDelta = Double.NaN;
        }
    }

    /**
     * Predicate to check that the disaster times are uniformly spaced within a
     * tolerance of {@link TimTam#timeEpsilon}. If there are less than two
     * disaster times then this defaults to false.
     */
    protected boolean disasterTimesAreUniform() {
        if (disasterTimesInput.get() == null) {
            return false;
        }

        Double[] dts = this.disasterTimesInput.get().getValues();
        if (dts.length < 2) {
            return false;
        } else if (dts.length == 2) {
            return true;
        } else {
            Double[] dds = new Double[dts.length-1];
            for (int ix=1; ix<dts.length; ix++) {
                dds[ix-1] = Math.abs(dts[ix-1] - dts[ix]);
            }
            List<Double> ddsL = Arrays.asList(dds);
            return (
                (Collections.max(ddsL) - Collections.min(ddsL)) < this.timeEpsilon
            );
        }
    }

    /**
     * Initialise the change time arrays.
     *
     * <p>Initialising the change time arrays is non-trivial because it depends
     * upon which parameterisation is being used. This method is not tuned for
     * performance so you may not want to call it within the MCMC loop.</p>
     *
     */
    private void initChangeTimes() {

        if (this.usingCanonical) {
            // This conditional checks that the number of change times is
            // appropriate for the number of lambda parameters that we are
            // estimating. This same pattern is repeated for the other parameters as
            // well.
            if (this.lambdaChangeTimesInput.get() != null) {
                checkChangeTimes(lambdaChangeTimesInput.get(),lambdaInput.get(),"lambda");
                this.lambdaChangeTimes = lambdaChangeTimesInput.get().getValues();
            } else {
                this.lambdaChangeTimes = new Double[]{};
            }

            if (this.muChangeTimesInput.get() != null) {
                checkChangeTimes(muChangeTimesInput.get(), muInput.get(), "mu");
                this.muChangeTimes = muChangeTimesInput.get().getValues();
            } else {
                this.muChangeTimes = new Double[]{};
            }

            if (this.psiChangeTimesInput.get() != null) {
                checkChangeTimes(psiChangeTimesInput.get(),psiInput.get(),"psi");
                this.psiChangeTimes = psiChangeTimesInput.get().getValues();
            } else {
                this.psiChangeTimes = new Double[]{};
            }

            if (this.omegaChangeTimesInput.get() != null) {
                checkChangeTimes(omegaChangeTimesInput.get(),omegaInput.get(),"omega");
                this.omegaChangeTimes = omegaChangeTimesInput.get().getValues();
            } else {
                this.omegaChangeTimes = new Double[]{};
            }

            if (rhoChangeTimesInput.get() != null) {
                this.rhoChangeTimes = rhoChangeTimesInput.get().getValues();
            } else {
                this.rhoChangeTimes = new Double[]{};
            }

            if (nuChangeTimesInput.get() != null) {
                this.nuChangeTimes = nuChangeTimesInput.get().getValues();
            } else {
                this.nuChangeTimes = new Double[]{};
            }
        } else if (this.usingR0) {

            int numLambdaChanges = 0;
            int numMuChanges = 0;
            int numPsiChanges = 0;
            int numOmegaChanges = 0;

            if (this.r0ChangeTimesInput.get() != null) {
                checkChangeTimes(r0ChangeTimesInput.get(),r0Input.get(),"r0");
                this.r0ChangeTimes = this.r0ChangeTimesInput.get().getValues();
                numLambdaChanges += this.r0ChangeTimes.length;
            } else {
                this.r0ChangeTimes = new Double[]{};
            }

            if (this.sigmaChangeTimesInput.get() != null) {
                checkChangeTimes(sigmaChangeTimesInput.get(),sigmaInput.get(),"sigma");
                this.sigmaChangeTimes = this.sigmaChangeTimesInput.get().getValues();
                numLambdaChanges += this.sigmaChangeTimes.length;
                numMuChanges += this.sigmaChangeTimes.length;
                numPsiChanges += this.sigmaChangeTimes.length;
                numOmegaChanges += this.sigmaChangeTimes.length;
            } else {
                this.sigmaChangeTimes = new Double[]{};
            }

            if (this.propPsiChangeTimesInput.get() != null) {
                checkChangeTimes(propPsiChangeTimesInput.get(),propPsiInput.get(),"p-sub-psi");
                this.propPsiChangeTimes = this.propPsiChangeTimesInput.get().getValues();
                numMuChanges += this.propPsiChangeTimes.length;
                numPsiChanges += this.propPsiChangeTimes.length;
            } else {
                this.propPsiChangeTimes = new Double[]{};
            }

            if (this.propOmegaChangeTimesInput.get() != null) {
                checkChangeTimes(propOmegaChangeTimesInput.get(),propOmegaInput.get(),"p-sub-omega");
                this.propOmegaChangeTimes = this.propOmegaChangeTimesInput.get().getValues();
                numMuChanges += this.propOmegaChangeTimes.length;
                numOmegaChanges += this.propOmegaChangeTimes.length;
            } else {
                this.propOmegaChangeTimes = new Double[]{};
            }

            if (numLambdaChanges == 0) {
                this.lambdaChangeTimes = new Double[]{};
            } else {
                this.lambdaChangeTimes = Numerics.concatenate(this.r0ChangeTimes,this.sigmaChangeTimes);
            }
            Arrays.sort(this.lambdaChangeTimes);

            if (numMuChanges == 0) {
                this.muChangeTimes = new Double[]{};
            } else {
                this.muChangeTimes = Numerics.concatenate(this.sigmaChangeTimes,this.propPsiChangeTimes,this.propOmegaChangeTimes);
            }
            Arrays.sort(this.muChangeTimes);

            if (numPsiChanges == 0) {
                this.psiChangeTimes = new Double[]{};
            } else {
                this.psiChangeTimes = Numerics.concatenate(this.sigmaChangeTimes,this.propPsiChangeTimes);
            }
            Arrays.sort(this.psiChangeTimes);

            if (numOmegaChanges == 0) {
                this.omegaChangeTimes = new Double[]{};
            } else {
                this.omegaChangeTimes = Numerics.concatenate(this.sigmaChangeTimes,this.propOmegaChangeTimes);
            }
            Arrays.sort(this.omegaChangeTimes);

            if (rhoChangeTimesInput.get() != null) {
                this.rhoChangeTimes = rhoChangeTimesInput.get().getValues();
            } else {
                this.rhoChangeTimes = new Double[]{};
            }

            if (nuChangeTimesInput.get() != null) {
                this.nuChangeTimes = nuChangeTimesInput.get().getValues();
            } else {
                this.nuChangeTimes = new Double[]{};
            }
        } else if (this.usingTimeSeriesR0) {
            int numLambdaChanges = 0;
            int numMuChanges = 0;
            int numPsiChanges = 0;
            int numNuChanges = 0;

            if (this.r0ChangeTimesInput.get() != null) {
                checkChangeTimes(r0ChangeTimesInput.get(),r0Input.get(),"r0");
                this.r0ChangeTimes = this.r0ChangeTimesInput.get().getValues();
                numLambdaChanges += this.r0ChangeTimes.length;
            } else {
                this.r0ChangeTimes = new Double[]{};
            }

            if (this.sigmaChangeTimesInput.get() != null) {
                checkChangeTimes(sigmaChangeTimesInput.get(),sigmaInput.get(),"sigma");
                this.sigmaChangeTimes = this.sigmaChangeTimesInput.get().getValues();
                numLambdaChanges += this.sigmaChangeTimes.length;
                numMuChanges += this.sigmaChangeTimes.length;
                numPsiChanges += this.sigmaChangeTimes.length;
                numNuChanges += this.sigmaChangeTimes.length;
            } else {
                this.sigmaChangeTimes = new Double[]{};
            }

            if (this.propPsiChangeTimesInput.get() != null) {
                checkChangeTimes(propPsiChangeTimesInput.get(),propPsiInput.get(),"p-sub-psi");
                this.propPsiChangeTimes = this.propPsiChangeTimesInput.get().getValues();
                numMuChanges += this.propPsiChangeTimes.length;
                numPsiChanges += this.propPsiChangeTimes.length;
            } else {
                this.propPsiChangeTimes = new Double[]{};
            }

            if (this.propTimeSeriesChangeTimesInput.get() != null) {
                checkChangeTimes(propTimeSeriesChangeTimesInput.get(),propTimeSeriesInput.get(),"p-sub-ts");
                this.propTimeSeriesChangeTimes = this.propTimeSeriesChangeTimesInput.get().getValues();
                numMuChanges += this.propTimeSeriesChangeTimes.length;
                numNuChanges += this.propTimeSeriesChangeTimes.length;
            } else {
                this.propTimeSeriesChangeTimes = new Double[]{};
            }

            if (numLambdaChanges == 0) {
                this.lambdaChangeTimes = new Double[]{};
            } else {
                // We need to be careful here, if one of the parameters does not
                // changes through time then the value will be null so we need
                // to use an empty array for its values. The temporary variables
                // have been used because then we can always just concatenate
                // the results, otherwise we would need to have a huge
                // conditional to handle every possible combination. The same
                // applies for the other parameters below.
                Double[] tmpR0ChangeTimes;
                Double[] tmpSigmaChangeTimes;
                if (this.r0ChangeTimesInput.get() != null) {
                    tmpR0ChangeTimes =
                        this.r0ChangeTimesInput.get().getValues();
                } else {
                    tmpR0ChangeTimes = new Double[]{};
                }
                if (this.sigmaChangeTimesInput.get() != null) {
                    tmpSigmaChangeTimes =
                        this.sigmaChangeTimesInput.get().getValues();
                } else {
                    tmpSigmaChangeTimes = new Double[]{};
                }
                this.lambdaChangeTimes = Numerics.concatenate(
                    tmpR0ChangeTimes,
                    tmpSigmaChangeTimes
                );
            }
            Arrays.sort(this.lambdaChangeTimes);

            if (numMuChanges == 0) {
                this.muChangeTimes = new Double[]{};
            } else {
                Double[] tmpSigmaChangeTimes =
                    this.sigmaChangeTimesInput.get() != null ?
                        this.sigmaChangeTimesInput.get().getValues() :
                        new Double[]{};
                Double[] tmpPropPsiChangeTimes =
                    this.propPsiChangeTimesInput.get() != null ?
                        this.propPsiChangeTimesInput.get().getValues() :
                        new Double[]{};
                Double[] tmpPropTSChangeTimes =
                    this.propTimeSeriesChangeTimesInput.get() != null ?
                        this.propTimeSeriesChangeTimesInput.get().getValues() :
                        new Double[]{};

                this.muChangeTimes = Numerics.concatenate(
                    tmpSigmaChangeTimes,
                    tmpPropPsiChangeTimes,
                    tmpPropTSChangeTimes
                );
            }
            Arrays.sort(this.muChangeTimes);

            if (numPsiChanges == 0) {
                this.psiChangeTimes = new Double[]{};
            } else {
                Double[] tmpSigmaChangeTimes =
                    this.sigmaChangeTimesInput.get() != null ?
                        this.sigmaChangeTimesInput.get().getValues() :
                        new Double[]{};
                Double[] tmpPropPsiChangeTimes =
                    this.propPsiChangeTimesInput.get() != null ?
                        this.propPsiChangeTimesInput.get().getValues() :
                        new Double[]{};
                this.psiChangeTimes = Numerics.concatenate(
                    tmpSigmaChangeTimes,
                    tmpPropPsiChangeTimes
                );
            }
            Arrays.sort(this.psiChangeTimes);

            if (numNuChanges == 0) {
                this.nuChangeTimes = new Double[]{};
            } else {
                Double[] tmpSigmaChangeTimes =
                    this.sigmaChangeTimesInput.get() != null ?
                        this.sigmaChangeTimesInput.get().getValues() :
                        new Double[]{};
                Double[] tmpPropTSChangeTimes =
                    this.propTimeSeriesChangeTimesInput.get() != null ?
                        this.propTimeSeriesChangeTimesInput.get().getValues() :
                        new Double[]{};
                this.nuChangeTimes = Numerics.concatenate(
                    tmpSigmaChangeTimes,
                    tmpPropTSChangeTimes
                );
            }
            Arrays.sort(this.nuChangeTimes);

            // Because this model assumes there is no rho or omega type
            // surveillance.
            this.rhoChangeTimes = new Double[]{};
            this.omegaChangeTimes = new Double[]{};
        }


        // Use a set to maintain a collection of the unique parameter change
        // times encountered to make it easier to handle the case where multiple
        // parameters change at the same point in time.
        Set<Double> uniqParamChangeTimes = new HashSet<>(Collections.emptySet());
        uniqParamChangeTimes.addAll(Arrays.asList(this.lambdaChangeTimes));
        uniqParamChangeTimes.addAll(Arrays.asList(this.muChangeTimes));
        uniqParamChangeTimes.addAll(Arrays.asList(this.psiChangeTimes));
        uniqParamChangeTimes.addAll(Arrays.asList(this.omegaChangeTimes));
        uniqParamChangeTimes.addAll(Arrays.asList(this.rhoChangeTimes));
        uniqParamChangeTimes.addAll(Arrays.asList(this.nuChangeTimes));
        if (uniqParamChangeTimes.size() > 0) {
            this.paramChangeTimes = uniqParamChangeTimes.toArray(new Double[]{0.0});
            Arrays.sort(this.paramChangeTimes);
        } else {
            this.paramChangeTimes = new Double[]{};
        }
    }

    /**
     * Throw an informative error message if something about the change times
     * given do not make sense.
     */
    private void checkChangeTimes(RealParameter changeTimes,
                                  RealParameter values,
                                  String name) {
        if (values.getDimension() != (changeTimes.getDimension() + 1)) {
            throw new RuntimeException(
                "The dimension of " + name + " and the number of times it changes are incompatible."
            );
        }
    }


    /**
     * Update the history sizes and times.
     *
     * <p>The history times should not change, so this is not entirely
     * necessary, unless there are uncertain tip dates on the tree in which case
     * this would all break.</p>
     */
    private boolean updateHistoryChecks() {
        boolean hasNegativeSize = false;
        if (this.numHistoryChecks > 0) {
            for (int ix = 0; ix < this.numHistoryChecks; ix++) {
                if (this.historySizesInput.get().getNativeValue(ix) >= 0) {
                    this.historySizes[ix] = this.historySizesInput.get().getNativeValue(ix);
                    this.historyTimes[ix] = this.historyTimesInput.get().getArrayValue(ix);
                } else {
                    hasNegativeSize = true;
                }
            }
        }
        return hasNegativeSize;
    }

    /**
     * <p>Predicate for being an unscheduled node.</p>
     *
     * <p>Because there is the possibility of numerical issues we use the
     * timeEpsilon to check if a node is sufficiently close to a catastrophe to
     * be considered one.</p>
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

    public double birth(int intNum) {
        return this.lambdaValues[intNum];
    }

    public double death(int intNum) {
        return this.muValues[intNum];
    }

    public double psi(int intNum) {
        return this.psiValues[intNum];
    }

    public double rho(int intNum) {
        return this.rhoValues[intNum];
    }

    public double omega(int intNum) {
        return this.omegaValues[intNum];
    }

    public double nu(int intNum) {
        return this.nuValues[intNum];
    }

    @Override
    public double calculateLogP() {
        // If the tree is tall enough that the root happens before the origin
        // then this point in parameter space has probability zero.
        if (this.tree.getRoot().getHeight() >= this.originTime.getValue()) {
            return Double.NEGATIVE_INFINITY;
        }

        this.nb.setIsDegenerate(0);

        // If any of the history sizes are negative, then the log-likelihood is
        // negative infinity. In this case we can avoid the rest of the
        // calculation and return this early to save time.
        if (updateHistoryChecks()) {
            return Double.NEGATIVE_INFINITY;
        }

        updateIntervalTerminators();
        updateLTTArray();
        updateRateAndProbParams();
        for (int ix = 0; ix < this.numTimeIntervals; ix++) {
            updateOdeHelpers(ix);
            processInterval(ix);
            processObservation(ix);
        }

        // if the likelihood conditions upon the observation of the process then
        // we need to account for this in the log-likelihood.
        if (this.conditionOnObservation) {
            if (this.paramChangeTimes.length == 0) {
            double probUnobserved = p0(this.timeFromOriginToFinalDatum, 0.0, 1.0);
            return Numerics.arraySum(this.lnls) + Numerics.arraySum(this.lncs) - Math.log(1 - probUnobserved);
            } else {
                throw new RuntimeException("conditioning upon observation is not yet implemented for varying rates.");
            }
        } else {
            return Numerics.arraySum(this.lnls) + Numerics.arraySum(this.lncs);
        }
    }

    /**
     * <p>This function updates the interval terminators which could change if
     * there are changes to the Tree object that lead to the order of events
     * changing.</p>
     *
     * <p>Note that all the times are relative to a zero which is the time of
     * the last sequenced sample, i.e. the most recent leaf in the tree. This
     * means that when parameter change times, and occurrence observations and
     * any other events are included, they need to be with relation to this
     * time.</p>
     *
     * @implNote The use of {@link java.util.Arrays#sort} here is appropriate according to the
     * <a href="https://docs.oracle.com/javase/7/docs/api/java/util/Arrays.html#sort(java.lang.Object[])">documentation</a>
     * because it is suitable for sorting concatenated sorted arrays.
     */
    private void updateIntervalTerminators() {
        int iTx = 0;

        for (Double paramChangeTime : this.paramChangeTimes) {
            this.intervalTerminators[iTx].setTypeTimeAndCount(
                    "paramValueChange",
                    paramChangeTime);
            iTx++;
        }

        for (Double occurrenceTime : this.occurrenceTimes) {
            this.intervalTerminators[iTx].setTypeTimeAndCount(
                    "occurrence",
                    occurrenceTime);
            iTx++;
        }

        for (int ix = 0; ix < this.catastropheTimes.length; ix++) {
            this.intervalTerminators[iTx].setTypeTimeAndCount(
                    "catastrophe",
                    this.catastropheTimes[ix],
                    this.catastropheSizes[ix]
            );
            iTx++;
        }

        for (int ix = 0; ix < this.disasterTimes.length; ix++) {
            this.intervalTerminators[iTx].setTypeTimeAndCount(
                    "disaster",
                    this.disasterTimes[ix],
                    this.disasterSizes[ix]
            );
            iTx++;
        }

        // Note that we get this values from the Tree object directly because
        // they may have changes since this method was last called.
        for (Node node : this.tree.getNodesAsArray()) {
            if (isUnscheduledTreeNode(node)) {
                this.intervalTerminators[iTx].setTypeTimeAndCount(
                        node.isLeaf() ? "sample" : "birth",
                        node.getHeight());
                iTx++;
            }
        }

        // Since this history checks array should have been updated prior to
        // calling the update function to update the interval terminators the
        // values being read here should be correct.
        for (int ix = 0; ix < this.historyTimes.length; ix++) {
            this.intervalTerminators[iTx].setTypeTimeAndCount(
                    "historyEstimate",
                    this.historyTimes[ix],
                    this.historySizes[ix]
            );
            iTx++;
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

    private void initIntervalTerminators() {
        this.intervalStartTimes = new double[this.numTimeIntervals];
        this.intervalEndTimes = new double[this.numTimeIntervals];

        this.intervalTerminators = new IntervalTerminator[this.numTimeIntervals];
        for (int ix = 0; ix < this.numTimeIntervals; ix++) {
            this.intervalTerminators[ix] = new IntervalTerminator();
        }
    }

    private void initRateAndProbParams() {
        this.lambdaValues = new double[this.numTimeIntervals];
        this.muValues = new double[this.numTimeIntervals];
        this.psiValues = new double[this.numTimeIntervals];
        this.rhoValues = new double[this.numTimeIntervals];
        this.omegaValues = new double[this.numTimeIntervals];
        this.nuValues = new double[this.numTimeIntervals];
    }

    private void initLTTArray() {
        this.kValues = new int[this.numTimeIntervals];
    }

    private void updateLTTArray() {
        String eventType;
        int crrK = 1;
        this.kValues[0] = crrK; // TODO This should be read from the XML!!!

        for (int ix=1; ix<this.numTimeIntervals; ix++) {
            eventType = this.intervalTerminators[ix-1].getType();
            if (eventType.equals("birth")) {
                crrK += 1;
            } else if (eventType.equals("sample")) {
                crrK -= 1;
            } else if (eventType.equals("catastrophe")) {
                crrK -= this.intervalTerminators[ix-1].getCount();
            } else {
                if (
                    !(
                        eventType.equals("occurrence") |
                            eventType.equals("paramValueChange") |
                            eventType.equals("disaster") |
                            eventType.equals("historyEstimate")
                    )
                ) {
                    throw new RuntimeException("Unexpected interval terminator type: " + eventType);
                }
            }
            this.kValues[ix] = crrK;
        }
    }

    /**
     * This function is called to update parameters of the model once one of the
     * inputs to the model has changed. Because it is called within the MCMC
     * loop it needs to be efficient.
     */
    private void updateRateAndProbParams() {
        Double nxtBwdTime;
        int ixLambda = 0;
        int ixMu = 0;
        int ixPsi = 0;
        int ixOmega = 0;
        int ixRho = 0;
        int ixNu = 0;
        Double crrLambda, crrMu, crrPsi, crrOmega, crrRho, crrNu;
        Double nxtLambdaChange,nxtMuChange,nxtPsiChange,
            nxtOmegaChange,nxtRhoChange,nxtNuChange;

        if (this.usingCanonical) {
            Double[] lambdas, mus, psis, omegas, rhos, nus;
            lambdas = this.lambdaInput.get().getValues();
            crrLambda = lambdas[ixLambda];
            if (lambdas.length > 1) {
                nxtLambdaChange = this.lambdaChangeTimes[ixLambda];
            } else {
                nxtLambdaChange = Double.NEGATIVE_INFINITY;
            }

            mus = this.muInput.get().getValues();
            crrMu = mus[ixMu];
            if (mus.length > 1) {
                nxtMuChange = this.muChangeTimes[ixMu];
            } else {
                nxtMuChange = Double.NEGATIVE_INFINITY;
            }

            psis = this.psiInput.get().getValues();
            crrPsi = psis[ixPsi];
            if (psis.length > 1) {
                nxtPsiChange = this.psiChangeTimes[ixPsi];
            } else {
                nxtPsiChange = Double.NEGATIVE_INFINITY;
            }

            // The omega, rho and nu values are handled very slightly
            // differently because they have different default values. After
            // these blocks they can be treated the same though in the loop
            // below.
            if (this.omegaInput.get() != null) {
                omegas = this.omegaInput.get().getValues();
                crrOmega = omegas[ixOmega];
                if (omegas.length > 1) {
                    nxtOmegaChange = this.omegaChangeTimes[ixOmega];
                } else {
                    nxtOmegaChange = Double.NEGATIVE_INFINITY;
                }
            } else {
                omegas = new Double[]{};
                crrOmega = 0.0;
                nxtOmegaChange = Double.NEGATIVE_INFINITY;
            }

            if (this.rhoInput.get() != null) {
                rhos = this.rhoInput.get().getValues();
                crrRho = rhos[ixRho];
                if (rhos.length > 1) {
                    nxtRhoChange = this.rhoChangeTimes[ixRho];
                } else {
                    nxtRhoChange = Double.NEGATIVE_INFINITY;
                }
            } else {
                rhos = new Double[]{};
                crrRho = 0.0;
                nxtRhoChange = Double.NEGATIVE_INFINITY;
            }

            if (this.nuInput.get() != null) {
                nus = this.nuInput.get().getValues();
                crrNu = nus[ixNu];
                if (nus.length > 1) {
                    nxtNuChange = this.nuChangeTimes[ixNu];
                } else {
                    nxtNuChange = Double.NEGATIVE_INFINITY;
                }
            } else {
                nus = new Double[]{};
                crrNu = 0.0;
                nxtNuChange = Double.NEGATIVE_INFINITY;
            }

            this.lambdaValues[0] = crrLambda;
            this.muValues[0] = crrMu;
            this.psiValues[0] = crrPsi;
            this.omegaValues[0] = crrOmega;
            this.rhoValues[0] = crrRho;
            this.nuValues[0] = crrNu;

            for (int ix=1; ix<this.numTimeIntervals; ix++) {
                nxtBwdTime = this.intervalStartTimes[ix];

                if (nxtBwdTime <= nxtLambdaChange) {
                    ixLambda++;
                    crrLambda = lambdas[ixLambda];
                    nxtLambdaChange =
                        ixLambda < this.lambdaChangeTimes.length ?
                            this.lambdaChangeTimes[ixLambda] :
                            Double.NEGATIVE_INFINITY;
                }
                if (nxtBwdTime <= nxtMuChange) {
                    ixMu++;
                    crrMu = mus[ixMu];
                    nxtMuChange =
                        ixMu < this.muChangeTimes.length ?
                            this.muChangeTimes[ixMu] :
                            Double.NEGATIVE_INFINITY;
                }
                if (nxtBwdTime <= nxtPsiChange) {
                    ixPsi++;
                    crrPsi = psis[ixPsi];
                    nxtPsiChange =
                        ixPsi < this.psiChangeTimes.length ?
                            this.psiChangeTimes[ixPsi] :
                            Double.NEGATIVE_INFINITY;
                }
                if (nxtBwdTime <= nxtOmegaChange) {
                    ixOmega++;
                    crrOmega = omegas[ixOmega];
                    nxtOmegaChange =
                        ixOmega < this.omegaChangeTimes.length ?
                            this.omegaChangeTimes[ixOmega] :
                            Double.NEGATIVE_INFINITY;
                }
                if (nxtBwdTime <= nxtRhoChange) {
                    ixRho++;
                    crrRho = rhos[ixRho];
                    nxtRhoChange =
                        ixRho < this.rhoChangeTimes.length ?
                            this.rhoChangeTimes[ixRho] :
                            Double.NEGATIVE_INFINITY;
                }
                if (nxtBwdTime <= nxtNuChange) {
                    ixNu++;
                    crrNu = nus[ixNu];
                    nxtNuChange =
                        ixNu < this.nuChangeTimes.length ?
                            this.nuChangeTimes[ixNu] :
                            Double.NEGATIVE_INFINITY;
                }

                this.lambdaValues[ix] = crrLambda;
                this.muValues[ix] = crrMu;
                this.psiValues[ix] = crrPsi;
                this.omegaValues[ix] = crrOmega;
                this.rhoValues[ix] = crrRho;
                this.nuValues[ix] = crrNu;
            }
        } else if (this.usingR0) {
            int ixR0 = 0;
            int ixSigma = 0;
            int ixPropPsi = 0;
            int ixPropOmega = 0;
            Double[] r0s, sigmas, propPsis, propOmegas, rhos, nus;
            Double crrR0, crrSigma, crrPropPsi, crrPropOmega;
            Double nxtR0Change, nxtSigmaChange, nxtPropPsiChange,
                nxtPropOmegaChange;

            r0s = this.r0Input.get().getValues();
            crrR0 = r0s[ixR0];
            nxtR0Change =
                r0s.length > 1 ?
                    this.r0ChangeTimes[ixR0] :
                    Double.NEGATIVE_INFINITY;

            sigmas = this.sigmaInput.get().getValues();
            crrSigma = sigmas[ixSigma];
            nxtSigmaChange =
                sigmas.length > 1 ?
                    this.sigmaChangeTimes[ixSigma] :
                    Double.NEGATIVE_INFINITY;

            propPsis = this.propPsiInput.get().getValues();
            crrPropPsi = propPsis[ixPropPsi];
            nxtPropPsiChange =
                propPsis.length > 1 ?
                    this.propPsiChangeTimes[ixPropPsi] :
                    Double.NEGATIVE_INFINITY;

            propOmegas =
                this.propOmegaInput.get() != null ?
                 this.propOmegaInput.get().getValues() :
                    new Double[]{0.0};
            crrPropOmega = propOmegas[ixPropOmega];
            nxtPropOmegaChange =
                propOmegas.length > 1 ?
                    this.propOmegaChangeTimes[ixPropOmega] :
                    Double.NEGATIVE_INFINITY;

            if (this.rhoInput.get() != null) {
                rhos = this.rhoInput.get().getValues();
                crrRho = rhos[ixRho];
                if (rhos.length > 1) {
                    nxtRhoChange = this.rhoChangeTimes[ixRho];
                } else {
                    nxtRhoChange = Double.NEGATIVE_INFINITY;
                }
            } else {
                rhos = new Double[]{};
                crrRho = 0.0;
                nxtRhoChange = Double.NEGATIVE_INFINITY;
            }

            if (this.nuInput.get() != null) {
                nus = this.nuInput.get().getValues();
                crrNu = nus[ixNu];
                if (nus.length > 1) {
                    nxtNuChange = this.nuChangeTimes[ixNu];
                } else {
                    nxtNuChange = Double.NEGATIVE_INFINITY;
                }
            } else {
                nus = new Double[]{};
                crrNu = 0.0;
                nxtNuChange = Double.NEGATIVE_INFINITY;
            }

            this.lambdaValues[0] = crrR0 * crrSigma;
            this.muValues[0] = (1 - crrPropPsi - crrPropOmega) * crrSigma;
            this.psiValues[0] = crrPropPsi * crrSigma;
            this.omegaValues[0] = crrPropOmega * crrSigma;
            this.rhoValues[0] = crrRho;
            this.nuValues[0] = crrNu;

            for (int ix=1; ix<this.numTimeIntervals; ix++) {
                nxtBwdTime = this.intervalStartTimes[ix];

                if (nxtBwdTime <= nxtR0Change) {
                    ixR0++;
                    crrR0 = r0s[ixR0];
                    nxtR0Change =
                        ixR0 < this.r0ChangeTimes.length ?
                            this.r0ChangeTimes[ixR0] :
                            Double.NEGATIVE_INFINITY;
                }

                if (nxtBwdTime <= nxtSigmaChange) {
                    ixSigma++;
                    crrSigma = sigmas[ixSigma];
                    nxtSigmaChange =
                        ixSigma < this.sigmaChangeTimes.length ?
                            this.sigmaChangeTimes[ixSigma] :
                            Double.NEGATIVE_INFINITY;
                }

                if (nxtBwdTime <= nxtPropPsiChange) {
                    ixPropPsi++;
                    crrPropPsi = propPsis[ixPropPsi];
                    nxtPropPsiChange =
                        ixPropPsi < this.propPsiChangeTimes.length ?
                            this.propPsiChangeTimes[ixPropPsi] :
                            Double.NEGATIVE_INFINITY;
                }

                if (nxtBwdTime <= nxtPropOmegaChange) {
                    ixPropOmega++;
                    crrPropOmega = propOmegas[ixPropOmega];
                    nxtPropOmegaChange =
                        ixPropOmega < this.propOmegaChangeTimes.length ?
                            this.propOmegaChangeTimes[ixPropOmega] :
                            Double.NEGATIVE_INFINITY;
                }

                if (nxtBwdTime <= nxtRhoChange) {
                    ixRho++;
                    crrRho = rhos[ixRho];
                    nxtRhoChange =
                        ixRho < this.rhoChangeTimes.length ?
                            this.rhoChangeTimes[ixRho] :
                            Double.NEGATIVE_INFINITY;
                }

                if (nxtBwdTime <= nxtNuChange) {
                    ixNu++;
                    crrNu = nus[ixNu];
                    nxtNuChange =
                        ixNu < this.nuChangeTimes.length ?
                            this.nuChangeTimes[ixNu] :
                            Double.NEGATIVE_INFINITY;
                }

                this.lambdaValues[ix] = crrR0 * crrSigma;
                this.muValues[ix] = (1 - crrPropPsi - crrPropOmega) * crrSigma;
                this.psiValues[ix] = crrPropPsi * crrSigma;
                this.omegaValues[ix] = crrPropOmega * crrSigma;
                this.rhoValues[ix] = crrRho;
                this.nuValues[ix] = crrNu;
            }
        } else if (this.usingTimeSeriesR0) {
            int ixR0 = 0;
            int ixSigma = 0;
            int ixPropPsi = 0;
            int ixPropTS = 0;

            Double[] r0s, sigmas, propPsis, omegas, propTSs, rhos, nus;
            Double crrR0, crrSigma, crrPropPsi, crrPropTS;
            Double nxtR0Change, nxtSigmaChange, nxtPropPsiChange,
                nxtPropTSChange;
            double tmp;

            r0s = this.r0Input.get().getValues();
            crrR0 = r0s[ixR0];
            nxtR0Change =
                r0s.length > 1 ?
                    this.r0ChangeTimes[ixR0] :
                    Double.NEGATIVE_INFINITY;

            sigmas = this.sigmaInput.get().getValues();
            crrSigma = sigmas[ixSigma];
            nxtSigmaChange =
                sigmas.length > 1 ?
                    this.sigmaChangeTimes[ixSigma] :
                    Double.NEGATIVE_INFINITY;

            propPsis = this.propPsiInput.get().getValues();
            crrPropPsi = propPsis[ixPropPsi];
            nxtPropPsiChange =
                propPsis.length > 1 ?
                    this.propPsiChangeTimes[ixPropPsi] :
                    Double.NEGATIVE_INFINITY;

            propTSs = this.propTimeSeriesInput.get().getValues();
            crrPropTS = propTSs[ixPropTS];
            nxtPropTSChange =
                propTSs.length > 1 ?
                    this.propTimeSeriesChangeTimes[ixPropTS] :
                    Double.NEGATIVE_INFINITY;

            if (this.rhoInput.get() != null | this.nuInput.get() != null) {
               throw new RuntimeException(
                   "The time series R0 parameterisation assumes rho and nu are null."
               ) ;
            } else {
                tmp = this.disasterTimeSeriesDelta * crrSigma * crrPropTS;
                crrNu = 2 * tmp / (2 + tmp);
            }

            this.lambdaValues[0] = crrR0 * crrSigma;
            this.muValues[0] = (1 - crrPropPsi - crrPropTS) * crrSigma;
            this.psiValues[0] = crrPropPsi * crrSigma;
            this.omegaValues[0] = 0.0;
            this.rhoValues[0] = 0.0;
            this.nuValues[0] = crrNu;

            for (int ix=1; ix<this.numTimeIntervals; ix++) {
                nxtBwdTime = this.intervalStartTimes[ix];

                if (nxtBwdTime <= nxtR0Change) {
                    ixR0++;
                    crrR0 = r0s[ixR0];
                    nxtR0Change =
                        ixR0 < this.r0ChangeTimes.length ?
                            this.r0ChangeTimes[ixR0] :
                            Double.NEGATIVE_INFINITY;
                }

                if (nxtBwdTime <= nxtSigmaChange) {
                    ixSigma++;
                    crrSigma = sigmas[ixSigma];
                    nxtSigmaChange =
                        ixSigma < this.sigmaChangeTimes.length ?
                            this.sigmaChangeTimes[ixSigma] :
                            Double.NEGATIVE_INFINITY;
                }

                if (nxtBwdTime <= nxtPropPsiChange) {
                    ixPropPsi++;
                    crrPropPsi = propPsis[ixPropPsi];
                    nxtPropPsiChange =
                        ixPropPsi < this.propPsiChangeTimes.length ?
                            this.propPsiChangeTimes[ixPropPsi] :
                            Double.NEGATIVE_INFINITY;
                }

                if (nxtBwdTime <= nxtPropTSChange) {
                    ixPropTS++;
                    crrPropTS = propTSs[ixPropTS];
                    nxtPropTSChange =
                        ixPropTS < this.propTimeSeriesChangeTimes.length ?
                            this.propTimeSeriesChangeTimes[ixPropTS] :
                            Double.NEGATIVE_INFINITY;
                }

                tmp = this.disasterTimeSeriesDelta * crrSigma * crrPropTS;
                crrNu = 2 * tmp / (2 + tmp);

                this.lambdaValues[ix] = crrR0 * crrSigma;
                this.muValues[ix] = (1 - crrPropPsi - crrPropTS) * crrSigma;
                this.psiValues[ix] = crrPropPsi * crrSigma;
                this.omegaValues[ix] = 0.0;
                this.rhoValues[ix] = 0.0;
                this.nuValues[ix] = crrNu;
            }
        }
    }

    /**
     * This method should mutate the attributes of this object to account for
     * the observation that occurred at the end of this interval.
     *
     * @implNote This method is implemented using if-else rather than a switch
     * because the switch syntax has changed between Java versions so this
     * guards against changes in versions.
     *
     * @param intNum the number of the interval ending in the observation.
     *
     * @see IntervalTerminator
     *
     */
    private void processObservation(int intNum) {
        IntervalTerminator intTerminator = this.intervalTerminators[intNum];
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
            // We need to store a value of zero here because otherwise, if the
            // order of observed events changes (for instance when the tree
            // changes) then values from previous calculations can sneak
            // through. Putting a zero here amounts to conditioning on the rate
            // to change at a fixed point in time.
            this.lnls[intNum] = 0.0;
        } else if (Objects.equals(intTypeStr, "historyEstimate")) {
            // We need to store the value of there actually being the estimated
            // number of hidden lineages at this point and then update the
            // distribution of hidden lineages to condition on there being that
            // many.
            int n = intTerminator.getCount();
            this.lnls[intNum] = this.nb.lnPMF(n);
            this.nb.setIsDegenerate(n);
        } else {
            throw new IllegalStateException(
                    "Unexpected value: " +
                    intTerminator.getType() +
                    "\n\tPlease see TimTamIntervalTerminator for types of intervals.");
        }
    }

    // this variable is just here in an attempt to resolve a memory leak...
    private final double[] tmpArry = {0,0,0,0,0};

    /**
     * Update the state of this object to account for the interval during which
     * there was no observation.
     *
     * <p>This method fills in an entry of the {@code this.lncs[intNum]} to
     * account for the duration of time that there was no observed event. See
     * Equation 1 of <a
     * href="https://doi.org/10.1371/journal.pcbi.1009805">Zarebski <emph>et
     * al</emph> (2022)</a> for further details. The method {@link
     * #processObservation(int)} is used to account for the actual observation
     * that was made at the end of the interval.</p>
     *
     * @param intNum the number of the interval to process.
     *
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

        // See Equation (12) from the Supporting information of the paper
        // mentioned above for rationale behind these expressions. The variables
        // represent the logarithm of the zero-th, first and second partial
        // derivatives of the PGF at the end of the interval. These are needed
        // to update the PGF so it represents the number of hidden lineages at
        // the end of the interval.
        double lnFM0;
        double lnFM1;
        double lnFM2;
        if (!this.nb.getIsZero()) {
            double lnPGFVal = this.nb.lnPGF(p0Val);
            double lnPGFDash1Val = this.nb.lnPGFDash1(p0Val);
            double lnPGFDash2Val = this.nb.lnPGFDash2(p0Val);
            if (k > 0) {
                lnFM0 = lnPGFVal + k * lnRVal;
                tmpArry[0] = lnPGFDash1Val + lnP0Dash1Val + k * lnRVal;
                tmpArry[1] = Math.log(k) + (k-1) * lnRVal + lnRDash1Val + lnPGFVal;
                lnFM1 = Numerics.logSumExp(tmpArry, 2);
                tmpArry[0] = lnPGFDash2Val + 2 * lnP0Dash1Val + k * lnRVal;
                tmpArry[1] = lnPGFDash1Val + lnP0Dash2Val + k * lnRVal;
                tmpArry[2] = Math.log(2) + lnPGFDash1Val + lnP0Dash1Val + Math.log(k) + (k - 1) * lnRVal + lnRDash1Val;
                tmpArry[3] = lnPGFVal + Math.log(k) + Math.log(k - 1) + (k - 2) * lnRVal + 2 * lnRDash1Val;
                tmpArry[4] = lnPGFVal + Math.log(k) + (k-1) * lnRVal + lnRDash2Val;
                lnFM2 = Numerics.logSumExp(tmpArry, 5);
            } else {
                lnFM0 = lnPGFVal;
                lnFM1 = lnPGFDash1Val + lnP0Dash1Val;
                tmpArry[0] = lnPGFDash2Val + 2 * lnP0Dash1Val;
                tmpArry[1] = lnPGFDash1Val + lnP0Dash2Val;
                lnFM2 = Numerics.logSumExp(tmpArry, 2);
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

    private double lnR(int intNum) {
        return Math.log(this.ohDiscriminant)
                + Math.log(this.ohExpFact)
                - (2 * Math.log(birth(intNum)))
                - (2 * Math.log((this.ohX2 - 1.0) - (this.ohX1 - 1.0) * this.ohExpFact));
    }

    private double lnRDash1(int intNum) {
        return Math.log(2)
                + Math.log(1 - this.ohExpFact)
                + Math.log(this.ohExpFact)
                + Math.log(this.ohDiscriminant)
                - 2 * Math.log(birth(intNum))
                - 3 * Math.log((this.ohX2 - this.ohExpFact * (this.ohX1 - 1.0) - 1.0));
    }

    private double lnRDash2(int intNum) {
        return Math.log(6)
                + 2 * Math.log(1 - this.ohExpFact)
                + Math.log(this.ohExpFact)
                + Math.log(this.ohDiscriminant)
                - 2 * Math.log(birth(intNum))
                - 4 * Math.log((this.ohX2 - this.ohExpFact * (this.ohX1 - 1.0) - 1.0));
    }

    protected double p0(double intervalDuration, double bwdTimeIntervalEnd, double z) {
        updateOdeHelpers(intervalDuration, bwdTimeIntervalEnd);
        return (this.ohX1 * (this.ohX2 - z) - this.ohX2 * (this.ohX1 - z) * this.ohExpFact) / ((this.ohX2 - z) - (this.ohX1 - z) * this.ohExpFact);
    }

    protected double p0(int intNum, double z) {
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

    protected double lnP0Dash1(int intNum) {
        double aa = this.ohX2 - this.ohX1 * this.ohExpFact;
        double bb = 1 - this.ohExpFact;
        double cc = this.ohX2 * this.ohExpFact - this.ohX1;
        return Math.log(cc * aa + this.ohX1 * this.ohX2 * Math.pow(bb, 2.0))
                - 2.0 * Math.log(aa - bb);
    }

    private double lnP0Dash2(int intNum) {
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
        int paramIx = 0;
        while (this.intervalEndTimes[paramIx] > bwdTimeIntervalEnd) {
            paramIx++;
        }
        double gamma = birth(paramIx) + death(paramIx) + psi(paramIx) + omega(paramIx);
        this.ohDiscriminant = Math.pow(gamma, 2.0) - 4.0 * birth(paramIx) * death(paramIx);
        double sqrtDisc = Math.sqrt(this.ohDiscriminant);
        this.ohX1 = (gamma - sqrtDisc) / (2 * birth(paramIx));
        this.ohX2 = (gamma + sqrtDisc) / (2 * birth(paramIx));
        this.ohExpFact = Math.exp(- sqrtDisc * intervalDuration);
    }

    /**
     * Update some expressions that are useful in multiple calculations. This
     * function gets called once before each interval is processed since the
     * results do not change within an interval.
     *
     * @implNote The local variables <code>br</code> and <code>dr</code> are
     * there to avoid repeated calls to the function that looks these values up.
     * The use of {@link java.lang.Math#pow} seems to be quicker than using
     * <code>gamma * gamma</code>, possibly because this is using native code.
     *
     * @param intNum the number of the interval to compute the ODE helper
     * functions for.
     */
    private void updateOdeHelpers(int intNum) {
        double intervalDuration = this.intervalStartTimes[intNum] - this.intervalEndTimes[intNum];
        double br = this.lambdaValues[intNum];
        double dr = this.muValues[intNum];
        double gamma = br + dr + psi(intNum) + omega(intNum);
        this.ohDiscriminant = Math.pow(gamma, 2.0) - 4.0 * br * dr;
        double sqrtDisc = Math.sqrt(this.ohDiscriminant);
        this.ohX1 = (gamma - sqrtDisc) / (2 * br);
        this.ohX2 = (gamma + sqrtDisc) / (2 * br);
        this.ohExpFact = Math.exp(- sqrtDisc * intervalDuration);
    }

    public HiddenLineageDist getTimTamNegBinom() {
        return this.nb;
    }


    @Override
    protected boolean requiresRecalculation() {
        return true;
    }

    /**
     * This is used to report the first interval duration as a way to assess if
     * the chain is mixing properly.
     *
     * @return the duration of the first interval.
     */
    protected double getFirstIntervalDuration() {
        return this.intervalStartTimes[0] - this.intervalEndTimes[0];
    }
}
