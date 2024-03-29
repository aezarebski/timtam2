package timtam;

import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.operator.ScaleOperator;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import org.junit.Test;

import java.util.function.BiPredicate;

import static org.junit.Assert.*;

public class TestTimTam {

    private final boolean verboseTesting = false;
    private final BiPredicate<Double, Double> approxEqual = (x, y) -> Math.abs(x - y) < 1e-5;
    private final BiPredicate<Double, Double> roughlyEqual = (x, y) -> Math.abs(x - y) < 1e-1;
    // check if within 5% of the second value.
    private final BiPredicate<Double, Double> kindaEqual = (x, y) -> (Math.abs(x - y) / Math.abs(y)) < 5e-2;

    /**
     * <p>This test draws on a similar one in BDSKY (see the
     * testLikelihoodCalculationSimple test in BirthDeathSkylineTest.java)
     * and checks that the TimTam values look similar in a special case.
     * These values were computed assuming constant parameters. The modified
     * version of testLikelihoodCalculationSimple in
     * BirthDeathSkylineTest.java is included at the end of this file.</p>
     *
     *  <table style="border:solid;">
     *  <tr>
     *  <th>R0</th>
     *  <th>lBDSKY</th>
     *  </tr>
     *  <tr>
     *  <td>1.5</td>
     *  <td>-26.10536013426608</td>
     *  </tr>
     *  <tr>
     *  <td>1.6</td>
     *  <td>-27.39912704449781</td>
     *  </tr>
     *  <tr>
     *  <td>1.7</td>
     *  <td>-28.76692080906782</td>
     *  </tr>
     *  <tr>
     *  <td>1.8</td>
     *  <td>-30.199269844913694</td>
     *  </tr>
     *  <tr>
     *  <td>1.9</td>
     *  <td>-31.68804261637826</td>
     *  </tr>
     *  <tr>
     *  <td>2.0</td>
     *  <td>-33.2262519906258</td>
     *  </tr>
     *  <tr>
     *  <td>3.0</td>
     *  <td>-50.33479549906616</td>
     *  </tr>
     *  <tr>
     *  <td>4.0</td>
     *  <td>-68.99855263104962</td>
     *  </tr>
     *  </table>
     */
    @Test
    public void testLikelihoodCalculationSimple() {

        TimTam tt = new TimTam();

        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);", false);
        tt.setInputValue("tree", tree);
        tt.setInputValue("originTime", new RealParameter("10.0"));

        double becomeUninfectiousRate = 1.5;
        double samplingProportion = 0.3;
        double[] r0Values = {1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 3.0, 4.0};
        double[] llhdValues = {
                -26.105360, -27.399127, -28.766920, -30.199269,
                -31.688042, -33.226251, -50.334795, -68.998552};

        tt.setInputValue("mu", new RealParameter(Double.toString(becomeUninfectiousRate * (1 - samplingProportion))));
        tt.setInputValue("psi", new RealParameter(Double.toString(becomeUninfectiousRate * samplingProportion)));

        tt.setInputValue("lambdaChangeTimes", null);
        String lambdaString;
        for (int ix = 0; ix < r0Values.length; ix++) {
            lambdaString = Double.toString(r0Values[ix] * becomeUninfectiousRate);
            tt.setInputValue("lambda", new RealParameter(lambdaString));
            tt.initAndValidate();
            assertTrue(kindaEqual.test(llhdValues[ix] - 2, tt.calculateLogP()));
        }
    }

    /**
     * <p>This test is similar the same as {@link TestTimTam#testLikelihoodCalculationSimple} except that it
     * includes usage of the {@link TimTam#historyTimesInput} to check
     * this is behaving as expected.</p>
     */
    @Test
    public void testLikelihoodCalculationSimpleWithHistory() {
        TimTam tt = new TimTam();
        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);", false);
        tt.setInputValue("tree", tree);
        tt.setInputValue("originTime", new RealParameter("10.0"));
        double becomeUninfectiousRate = 1.5;
        double samplingProportion = 0.3;
        tt.setInputValue("mu", new RealParameter(Double.toString(becomeUninfectiousRate * (1 - samplingProportion))));
        tt.setInputValue("psi", new RealParameter(Double.toString(becomeUninfectiousRate * samplingProportion)));

        double[] r0Values = {1.01, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 3.0, 4.0, 10.0};

        tt.setInputValue("lambdaChangeTimes", null);
        String lambdaString;

        int numSizes = 100;
        int[] historySizes = new int[numSizes];
        for (int i = 0; i < historySizes.length; i++) {
            historySizes[i] = i + 1;
        }
        double[] partialLlhdValues = new double[numSizes];
        double tmpLlhd;
        for (int ix = 0; ix < r0Values.length; ix++) {
            tt.setInputValue("historyTimes", null);
            tt.setInputValue("historySizes", null);
            lambdaString = Double.toString(r0Values[ix] * becomeUninfectiousRate);
            tt.setInputValue("lambda", new RealParameter(lambdaString));
            tt.initAndValidate();
            tmpLlhd = tt.calculateLogP();

            // Loop over the possible history sizes to see that
            // marginalising over them gives the same result.
            tt.setInputValue("historyTimes", new RealParameter("1.0"));
            for (int jx = 0; jx < historySizes.length; jx++) {
                tt.setInputValue("historySizes", new IntegerParameter(Integer.toString(historySizes[jx])));
                tt.initAndValidate();
                partialLlhdValues[jx] = tt.calculateLogP();
            }
            if (verboseTesting) {
                System.out.println("lambda = " + lambdaString);
                System.out.println(tmpLlhd);
                System.out.println(Numerics.logSumExp(partialLlhdValues));
            }
            assertTrue(kindaEqual.test(tmpLlhd, Numerics.logSumExp(partialLlhdValues)));

            tt.setInputValue("historySizes", new IntegerParameter("-1"));
            tt.initAndValidate();
            if (verboseTesting) {
                System.out.println("Now testing negative history size...");
            }
            assertTrue(tt.calculateLogP() == Double.NEGATIVE_INFINITY);
        }
    }

    /**
     * Helper function to convert numbers to parameters that BEAST likes.
     */
    private RealParameter asRealParam(double d) {
        return new RealParameter(Double.toString(d));
    }

    private RealParameter asRealParam(String s) {
        return new RealParameter(s);
    }

    /**
     * This test is a copy of {@link TestTimTam#testLikelihoodCalculationSimple}
     * where instead of using the canonical parameterisation we will use the r0
     * parameterisation.
     */
    @Test
    public void testLikelihoodCalculationSimpleR0() {

        double eps = 0.4;

        TimTam tt = new TimTam();

        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);", false);
        tt.setInputValue("tree", tree);
        tt.setInputValue("originTime", asRealParam("10.0"));

        double netBecomeUninfectiousRate = 1.5;
        double propSeq = 0.3;
        double propOcc = 0.0;
        double[] r0Values = {1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 3.0, 4.0};
        double[] llhdValues = {
                -26.105360, -27.399127, -28.766920, -30.199269,
                -31.688042, -33.226251, -50.334795, -68.998552};
        double[] canonicalParamValues = new double[llhdValues.length];

        for (int ix = 0; ix < r0Values.length; ix++) {
            tt.setInputValue("mu", asRealParam(netBecomeUninfectiousRate * (1 - propSeq - propOcc)));
            tt.setInputValue("psi", asRealParam(netBecomeUninfectiousRate * propSeq));
            tt.setInputValue("omega", asRealParam(netBecomeUninfectiousRate * propOcc));
            tt.setInputValue("lambda", asRealParam(r0Values[ix] * netBecomeUninfectiousRate));
            tt.setInputValue("parameterisation", "canonical");
            tt.initAndValidate();
            canonicalParamValues[ix] = tt.calculateLogP();
            assertEquals(llhdValues[ix] - 2, tt.calculateLogP(), eps);
        }

        tt.setInputValue("mu", null);
        tt.setInputValue("psi", null);
        tt.setInputValue("omega", null);
        tt.setInputValue("lambda", null);
        for (int ix = 0; ix < r0Values.length; ix++) {
            // basic reproduction number
            tt.setInputValue("r0",
                    new RealParameter(Double.toString(r0Values[ix])));
            // net removal rate
            tt.setInputValue("sigma",
                    new RealParameter(Double.toString(netBecomeUninfectiousRate)));
            // proportion observed and sequenced
            tt.setInputValue("propPsi",
                    new RealParameter(Double.toString(propSeq)));
            // proportion observed and not sequenced
            tt.setInputValue("propOmega",
                    new RealParameter(Double.toString(propOcc)));
            tt.setInputValue("parameterisation",
                    "r0");
            tt.initAndValidate();
            assertEquals(llhdValues[ix] - 2, tt.calculateLogP(), eps);
            assertEquals(canonicalParamValues[ix], tt.calculateLogP(), 1e-8);

        }
    }

    /**
     * <p>This test ensures that if an operator drops the origin below the height of the tree then the resulting log
     * probability becomes negative infinity.</p>
     */
    @Test
    public void testOriginLessThanHeightIsImpossible() {

        TimTam tt = new TimTam();
        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);", false);

        tt.setInputValue("tree", tree);
        RealParameter oT = new RealParameter("10.0");
        oT.setUpper(11.0);
        oT.setLower(0.1);
        tt.setInputValue("originTime", oT);

        tt.setInputValue("lambda", "2.25");
        tt.setInputValue("mu", "1.05");
        tt.setInputValue("psi", "0.45");

        tt.initAndValidate();

        assertTrue(kindaEqual.test(-26.105360134266082 - 2, tt.calculateLogP()));

        ScaleOperator operator = new ScaleOperator();
        operator.initByName(
                "parameter", oT,
                "scaleFactor", 0.5,
                "scaleAllIndependently", true,
                "weight", 10.0
        );

        // while the origin is older than the root of the tree the log-probability should be finite but once it drops
        // below that it should be infinite.
        while (oT.getValue() > 5.0) {
            assertTrue(tt.calculateLogP() > Double.NEGATIVE_INFINITY);
            operator.proposal();
        }
        assertTrue(tt.calculateLogP() == Double.NEGATIVE_INFINITY);
    }

    /**
     * This is pretty much the same as {@link #testLikelihoodCalculationSimple()} but
     * it uses a variable birth rate. The resulting log-likelihood should sit somewhere
     * between the likelihoods of the two extremes used.
     */
    @Test
    public void testLikelihoodCalculationVariableBirthRate() {

        TimTam tt = new TimTam();

        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);", false);
        tt.setInputValue("tree", tree);
        tt.setInputValue("originTime", new RealParameter("10.0"));
        tt.setInputValue("conditionOnObservation", false);

        double becomeUninfectiousRate = 1.5;
        double samplingProportion = 0.3;

        tt.setInputValue("mu", new RealParameter(Double.toString(becomeUninfectiousRate * (1 - samplingProportion))));
        tt.setInputValue("psi", new RealParameter(Double.toString(becomeUninfectiousRate * samplingProportion)));

        // We can check that when using a variable rate it interpolates between the extremes used.
        String lambdaString1p8 = Double.toString(1.8 * becomeUninfectiousRate);
        tt.setInputValue("lambda", new RealParameter(lambdaString1p8));
        tt.initAndValidate();
        double tmpWith1p8 = tt.calculateLogP();
        String lambdaString1p9 = Double.toString(1.9 * becomeUninfectiousRate);
        tt.setInputValue("lambda", new RealParameter(lambdaString1p9));
        tt.initAndValidate();
        double tmpWith1p9 = tt.calculateLogP();
        String lambdaStringV = lambdaString1p8 + " " + lambdaString1p9;
        tt.setInputValue("lambda", new RealParameter(lambdaStringV));
        tt.setInputValue("lambdaChangeTimes", new RealParameter("3.0"));
        tt.initAndValidate();
        double tmpWithVarying = tt.calculateLogP();
        assertTrue(tmpWith1p8 > tmpWithVarying);
        assertTrue(tmpWithVarying > tmpWith1p9);
    }

    /**
     * This is pretty much the same as {@link #testLikelihoodCalculationSimple()} but
     * it uses a variable death rate. The resulting log-likelihood should sit somewhere
     * between the likelihoods of the two extremes used. This is a little more complicated than
     * {@link #testLikelihoodCalculationVariableBirthRate()} because the becoming
     * uninfectious rate moves through the calculation.
     */
    @Test
    public void testLikelihoodCalculationVariableDeathRate() {

        TimTam tt = new TimTam();

        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);", false);
        tt.setInputValue("tree", tree);
        tt.setInputValue("originTime", new RealParameter("10.0"));
        tt.setInputValue("conditionOnObservation", false);
        tt.setInputValue("psi", new RealParameter("0.5"));
        tt.setInputValue("lambda", new RealParameter("2.5"));

        double lnPWith1p5, lnPWith1p6, lnPWithVarying;

        tt.setInputValue("mu", new RealParameter("1.5"));
        tt.initAndValidate();
        lnPWith1p5 = tt.calculateLogP();

        tt.setInputValue("mu", new RealParameter("1.6"));
        tt.initAndValidate();
        lnPWith1p6 = tt.calculateLogP();

        tt.setInputValue("mu", new RealParameter("1.5 1.6"));
        tt.setInputValue("muChangeTimes", new RealParameter("3.0"));
        tt.initAndValidate();
        lnPWithVarying = tt.calculateLogP();

        assertTrue((
                (lnPWith1p5 < lnPWithVarying & lnPWithVarying < lnPWith1p6) |
                        (lnPWith1p5 > lnPWithVarying & lnPWithVarying > lnPWith1p6)
        ));
    }

    @Test
    public void testLikelihood() {

        Tree tree = new TreeParser("(((1:3,2:1):1,3:4):2,4:6);", false);

        RealParameter catastropheTimes = new RealParameter();
        catastropheTimes.initByName("value", "0.0");

        RealParameter points = new RealParameter();
        points.initByName("value", "5.0 1.0");

        TimTam tt = new TimTam();
        tt.setInputValue("lambda", "2.0");
        tt.setInputValue("mu", "1.0");
        tt.setInputValue("psi", "0.3");
        tt.setInputValue("rho", "0.5");
        tt.setInputValue("omega", "0.6");
        tt.setInputValue("originTime", "7.0");
        tt.setInputValue("catastropheTimes", catastropheTimes);
        tt.setInputValue("occurrenceTimes", points);
        tt.setInputValue("tree", tree);
        tt.setInputValue("conditionOnObservation", "false");
        tt.initAndValidate();


        double fx, fxh, h, fxDash;
        h = 1e-6;
        fxh = tt.p0(1.0, 0.0, 1.0 - h);
        fx = tt.p0(1.0, 0.0, 1.0);
        fxDash = tt.lnP0Dash1(1.0, 0.0);
        assertTrue(approxEqual.test(Math.log((fx - fxh) / h), fxDash));

        assertTrue(
                roughlyEqual.test(-47.0, tt.calculateLogP()));

        double lnMean0 = tt.getTimTamNegBinom().getLnMean();
        tt.setInputValue("lambda", "3.0");
        tt.initAndValidate();
        assertFalse(
                roughlyEqual.test(-47.0, tt.calculateLogP())
        );
        double lnMean1 = tt.getTimTamNegBinom().getLnMean();
        assertFalse(roughlyEqual.test(lnMean0, lnMean1));

        // if we repeat this using an instance that conditions upon observation then the value should be different.
        TimTam ttConditioned = new TimTam();
        ttConditioned.setInputValue("lambda", "2.0");
        ttConditioned.setInputValue("mu", "1.0");
        ttConditioned.setInputValue("psi", "0.3");
        ttConditioned.setInputValue("rho", "0.5");
        ttConditioned.setInputValue("omega", "0.6");
//        ttConditioned.setInputValue("rootLength", "1.0");
        ttConditioned.setInputValue("originTime", "7.0");
        ttConditioned.setInputValue("catastropheTimes", catastropheTimes);
//        ttConditioned.setInputValue("points", points);
        ttConditioned.setInputValue("occurrenceTimes", points);
        ttConditioned.setInputValue("tree", tree);
//        ttConditioned.setInputValue("conditionOnObservation", "true");
        ttConditioned.initAndValidate();
        assertFalse(
                roughlyEqual.test(
                        -47.0,
                        ttConditioned.calculateLogP()));

    }

    @Test
    public void testVariableRhoProb() {
        String newickString =
                "(((1:2.0,2:1.0):2.0,3:4.0):2.0,((4:2.0,5:3.0):2.0,6:4.0):1.0);";
        Tree tree = new TreeParser(newickString, false);
        TimTam tt = new TimTam();
        tt.setInputValue("lambda", "2.0");
        tt.setInputValue("mu", "1.0");
        tt.setInputValue("psi", "0.0");
        tt.setInputValue("originTime", "7.0");
        tt.setInputValue("tree", tree);
        tt.setInputValue("conditionOnObservation", "false");
        tt.setInputValue("catastropheTimes", "1.0 0.0");

        double[] y = new double[4];

        tt.setInputValue("rho", "0.3 0.3");
        tt.setInputValue("rhoChangeTimes", "0.5");
        tt.initAndValidate();
        y[0] = tt.calculateLogP();

        tt.setInputValue("rho", "0.3 0.4");
        tt.setInputValue("rhoChangeTimes", "0.5");
        tt.initAndValidate();
        y[1] = tt.calculateLogP();

        tt.setInputValue("rho", "0.4 0.4");
        tt.setInputValue("rhoChangeTimes", "0.5");
        tt.initAndValidate();
        y[2] = tt.calculateLogP();

        tt.setInputValue("rho", "0.4");
        tt.initAndValidate();
        y[3] = tt.calculateLogP();

        assertTrue(y[0] != y[1]);
        assertTrue(
                (y[0] < y[1] & y[1] < y[2]) |
                        (y[2] < y[1] & y[1] < y[0])
        );
        assertTrue(approxEqual.test(y[2], y[3]));
    }

    @Test
    public void testVariableNuProb() {
        String newickString =
                "((1:3.0,2:1.0):3.0,3:1.0);";
        Tree tree = new TreeParser(newickString, false);
        TimTam tt = new TimTam();
        tt.setInputValue("lambda", "2.0");
        tt.setInputValue("mu", "1.0");
        tt.setInputValue("psi", "1.0");
        tt.setInputValue("originTime", "7.0");
        tt.setInputValue("tree", tree);
        tt.setInputValue("conditionOnObservation", "false");

        tt.setInputValue("disasterTimes", "4.0 1.0");
        tt.setInputValue("disasterSizes", "2 3");

        double[] y = new double[4];

        tt.setInputValue("nu", "0.3 0.3");
        tt.setInputValue("nuChangeTimes", "2.5");
        tt.initAndValidate();
        y[0] = tt.calculateLogP();

        tt.setInputValue("nu", "0.3 0.4");
        tt.setInputValue("nuChangeTimes", "2.5");
        tt.initAndValidate();
        y[1] = tt.calculateLogP();

        tt.setInputValue("nu", "0.4 0.4");
        tt.setInputValue("nuChangeTimes", "2.5");
        tt.initAndValidate();
        y[2] = tt.calculateLogP();

        tt.setInputValue("nu", "0.4");
        tt.initAndValidate();
        y[3] = tt.calculateLogP();

        assertTrue(y[0] != y[1]);
        assertTrue(
                (y[0] < y[1] & y[1] < y[2]) |
                        (y[2] < y[1] & y[1] < y[0])
        );
        assertTrue(approxEqual.test(y[2], y[3]));
    }

    @Test(expected = RuntimeException.class)
    public void testTMRCABeforeOriginThrowsException() {

        Tree tree = new TreeParser("((1:6.0,2:4.0):2.0,3:4.0);", false);
        TimTam tt = new TimTam();
        tt.setInputValue("lambda", "3.0");
        tt.setInputValue("mu", "1.0");
        tt.setInputValue("psi", "1.0");
        tt.setInputValue("originTime", "5.0");
        tt.setInputValue("tree", tree);
        tt.setInputValue("conditionOnObservation", "false");
        tt.initAndValidate();
    }

    @Test(expected = RuntimeException.class)
    public void testInvalidDimensionThrowsExceptionA() {

        Tree tree = new TreeParser("((1:6.0,2:4.0):2.0,3:4.0);", false);
        TimTam tt = new TimTam();
        tt.setInputValue("lambda", "3.0 4.0 5.0");
        tt.setInputValue("lambdaChangeTimes", "1.0");
        tt.setInputValue("mu", "1.0");
        tt.setInputValue("psi", "1.0");
        tt.setInputValue("originTime", "10.0");
        tt.setInputValue("tree", tree);
        tt.setInputValue("conditionOnObservation", "false");
        tt.initAndValidate();
    }

    @Test(expected = RuntimeException.class)
    public void testInvalidDimensionThrowsExceptionB() {

        Tree tree = new TreeParser("((1:6.0,2:4.0):2.0,3:4.0);", false);
        TimTam tt = new TimTam();
        tt.setInputValue("lambda", "3.0");
        tt.setInputValue("lambdaChangeTimes", "1.0");
        tt.setInputValue("mu", "1.0");
        tt.setInputValue("psi", "1.0");
        tt.setInputValue("originTime", "10.0");
        tt.setInputValue("tree", tree);
        tt.setInputValue("conditionOnObservation", "false");
        tt.initAndValidate();
    }

    @Test
    public void testDisasterTimesAreUniform() {
        TimTam tt = new TimTam();

        assertFalse(tt.disasterTimesAreUniform());

        tt.setInputValue("disasterTimes", "4.0");
        assertFalse(tt.disasterTimesAreUniform());

        tt.setInputValue("disasterTimes", "4.0 1.0");
        assertTrue(tt.disasterTimesAreUniform());

        tt.setInputValue("disasterTimes", "7.0 4.0 1.0");
        assertTrue(tt.disasterTimesAreUniform());

        tt.setInputValue("disasterTimes",
                Double.toString(7.0 + 0.1 * tt.timeEpsilon) + " 4.0 1.0");
        assertTrue(tt.disasterTimesAreUniform());

        tt.setInputValue("disasterTimes", "7.01 4.0 1.0");
        assertFalse(tt.disasterTimesAreUniform());

        tt.setInputValue("disasterTimes", "6.99 4.0 1.0");
        assertFalse(tt.disasterTimesAreUniform());
    }
}


// public void testLikelihoodCalculationSimple() throws Exception {

//     BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

//     Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
//     bdssm.setInputValue("tree", tree);
//     bdssm.setInputValue("origin", new RealParameter("10."));
//     bdssm.setInputValue("conditionOnSurvival", false);
//     bdssm.setInputValue("removalProbability", "1");

//     bdssm.setInputValue("reproductiveNumber", new RealParameter("1.5"));
//     bdssm.setInputValue("becomeUninfectiousRate", new RealParameter("1.5"));
//     bdssm.setInputValue("samplingProportion", new RealParameter("0.3") );

//     bdssm.initAndValidate();
//     bdssm.printTempResults = false;

//     assertEquals(-26.105360134266082, bdssm.calculateTreeLogLikelihood(tree), 1e-4);

//     System.out.println("| R0  | lBDSKY             |");
//     System.out.println("|-----+--------------------|");
//     String[] reproductiveNumbers =
//         {"1.5", "1.6", "1.7", "1.8", "1.9", "2.0", "3.0", "4.0"};
//     for (String r0String : reproductiveNumbers) {
//         bdssm.setInputValue("reproductiveNumber", new RealParameter(r0String));
//         bdssm.initAndValidate();
//         double tmp = bdssm.calculateTreeLogLikelihood(tree);
//         System.out.println("| " + r0String + " | " + tmp + " |");
//     }
// }
