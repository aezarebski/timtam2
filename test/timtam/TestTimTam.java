package timtam;

import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import org.junit.Test;

import java.util.function.BiPredicate;

import static org.junit.Assert.*;

public class TestTimTam {

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

        TimTam tt =  new TimTam();

        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
        tt.setInputValue("tree", tree);
        tt.setInputValue("originTime", new RealParameter("10.0"));

        double becomeUninfectiousRate = 1.5;
        double samplingProportion = 0.3;
        double[] r0Values = {1.5,1.6,1.7,1.8,1.9,2.0,3.0,4.0};
        double[] llhdValues = {
                -26.105360,-27.399127,-28.766920,-30.199269,
                -31.688042,-33.226251,-50.334795,-68.998552};

        tt.setInputValue("mu",new RealParameter(Double.toString(becomeUninfectiousRate * (1 - samplingProportion))));
        tt.setInputValue("psi",new RealParameter(Double.toString(becomeUninfectiousRate * samplingProportion)));

        tt.setInputValue("lambdaChangeTimes", null);
        String lambdaString;
        for (int ix = 0; ix < r0Values.length; ix++) {
            lambdaString = Double.toString(r0Values[ix] * becomeUninfectiousRate);
            tt.setInputValue("lambda", new RealParameter(lambdaString));
            tt.initAndValidate();
            assertTrue(kindaEqual.test(llhdValues[ix] - 2, tt.calculateLogP()));
        }
    }

    /** This is pretty much the same as {@link #testLikelihoodCalculationSimple()} but
     * it uses a variable birth rate. The resulting log-likelihood should sit somewhere
     * between the likelihoods of the two extremes used.
     */
    @Test
    public void testLikelihoodCalculationVariableBirthRate() {

        TimTam tt =  new TimTam();

        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
        tt.setInputValue("tree", tree);
        tt.setInputValue("originTime", new RealParameter("10.0"));
        tt.setInputValue("conditionOnObservation", false);

        double becomeUninfectiousRate = 1.5;
        double samplingProportion = 0.3;

        tt.setInputValue("mu",new RealParameter(Double.toString(becomeUninfectiousRate * (1 - samplingProportion))));
        tt.setInputValue("psi",new RealParameter(Double.toString(becomeUninfectiousRate * samplingProportion)));

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

    /** This is pretty much the same as {@link #testLikelihoodCalculationSimple()} but
     * it uses a variable death rate. The resulting log-likelihood should sit somewhere
     * between the likelihoods of the two extremes used. This is a little more complicated than
     * {@link #testLikelihoodCalculationVariableBirthRate()} because the becoming
     * uninfectious rate moves through the calculation.
     */
    @Test
    public void testLikelihoodCalculationVariableDeathRate() {

        TimTam tt =  new TimTam();

        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
        tt.setInputValue("tree", tree);
        tt.setInputValue("originTime", new RealParameter("10.0"));
        tt.setInputValue("conditionOnObservation", false);
        tt.setInputValue("psi",new RealParameter("0.5"));
        tt.setInputValue("lambda",new RealParameter("2.5"));

        double lnPWith1p5, lnPWith1p6, lnPWithVarying;

        tt.setInputValue("mu",new RealParameter("1.5"));
        tt.initAndValidate();
        lnPWith1p5 = tt.calculateLogP();

        tt.setInputValue("mu",new RealParameter("1.6"));
        tt.initAndValidate();
        lnPWith1p6 = tt.calculateLogP();

        tt.setInputValue("mu",new RealParameter("1.5 1.6"));
        tt.setInputValue("muChangeTimes",new RealParameter("3.0"));
        tt.initAndValidate();
        lnPWithVarying = tt.calculateLogP();

        assertTrue((
                (lnPWith1p5 < lnPWithVarying & lnPWithVarying < lnPWith1p6) |
                (lnPWith1p5 > lnPWithVarying & lnPWithVarying > lnPWith1p6)
        ));
    }

    @Test
    public void testLikelihood() {

        Tree tree = new TreeParser("(((1:3,2:1):1,3:4):2,4:6);",false);

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
                roughlyEqual.test(-47.0,tt.calculateLogP()));

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
    public void testLogSumExp() {

        // > log(sum(exp(c(1.17,0.95,0.10,1.58))))
        // [1] 2.465369
        // > log(sum(exp(c(1.17,0.95,0.10))))
        // [1] 1.933385

        TimTam tt = new TimTam();

        double[] xs = {1.17,0.95,0.10,1.58};

        assertTrue(approxEqual.test(2.465369, tt.logSumExp(xs, 4)));
        assertTrue(approxEqual.test(1.933385, tt.logSumExp(xs, 3)));


        // > log(sum(exp(c(-1.0,0.95,0.10,1.58))))
        // [1] 2.187591
        xs[0] = -1.0;
        assertTrue(approxEqual.test(2.187591, tt.logSumExp(xs)));
        // > log(sum(exp(c(1.0,0.95,0.10,1.58))))
        // [1] 2.421622
        xs[0] = 1.0;
        assertTrue(approxEqual.test(2.421622, tt.logSumExp(xs)));
        // > log(sum(exp(c(10.0,0.95,0.10,1.58))))
        // [1] 10.00039
        xs[0] = 10.0;
        assertTrue(approxEqual.test(10.00039, tt.logSumExp(xs)));
    }

    @Test
    public void testVariableBirthRate() {

        Tree tree = new TreeParser("(((1:3,2:1):1,3:4):2,4:6);",false);
        RealParameter occTimes = new RealParameter();
        occTimes.initByName("value", "5.0 1.0");

        double lv0 = 3.1;
        double lv1 = 4.1;
        double lv2 = 5.9;
        double lct0 = 4.3;
        double lct1 = 2.1;
        TimTam tt = new TimTam();
        tt.setInputValue("lambda", lv0 + " " + lv1 + " " + lv2);
        tt.setInputValue("lambdaChangeTimes", lct0 + " " + lct1);
        tt.setInputValue("mu", "1.0");
        tt.setInputValue("psi", "0.3");
        tt.setInputValue("rho", "0.5");
        tt.setInputValue("omega", "0.6");
        tt.setInputValue("originTime", "7.0");
        tt.setInputValue("occurrenceTimes", occTimes);
        tt.setInputValue("tree", tree);
        tt.setInputValue("conditionOnObservation", "false");
        tt.initAndValidate();

        assertTrue(approxEqual.test(tt.birth(lct0 + 0.1), lv0));
        assertTrue(approxEqual.test(tt.birth(lct0), lv1));
        assertTrue(approxEqual.test(tt.birth(lct0 - 0.1), lv1));

        assertTrue(approxEqual.test(tt.birth(lct1 + 0.1), lv1));
        assertTrue(approxEqual.test(tt.birth(lct1), lv2));
        assertTrue(approxEqual.test(tt.birth(lct1 - 0.1), lv2));
    }

    @Test
    public void testVariableRhoProb() {
        String newickString =
            "(((1:2.0,2:1.0):2.0,3:4.0):2.0,((4:2.0,5:3.0):2.0,6:4.0):1.0);";
        Tree tree = new TreeParser(newickString,false);
        TimTam tt = new TimTam();
        tt.setInputValue("lambda", "2.0");
        tt.setInputValue("mu", "1.0");
        tt.setInputValue("psi", "0.0");
        tt.setInputValue("originTime", "7.0");
        tt.setInputValue("tree", tree);
        tt.setInputValue("conditionOnObservation", "false");

        double[] y = new double[3];

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

        assertTrue(y[0] != y[1]);
        assertTrue(
                (y[0] < y[1] & y[1] < y[2]) |
                        (y[2] < y[1] & y[1] < y[0])
        );
    }

    @Test(expected = RuntimeException.class)
    public void testTMRCABeforeOriginThrowsException() {

        Tree tree = new TreeParser("((1:6.0,2:4.0):2.0,3:4.0);",false);
        TimTam tt = new TimTam();
        tt.setInputValue("lambda", "3.0");
        tt.setInputValue("mu", "1.0");
        tt.setInputValue("psi", "1.0");
        tt.setInputValue("originTime", "5.0");
        tt.setInputValue("tree", tree);
        tt.setInputValue("conditionOnObservation", "false");
        tt.initAndValidate();
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
