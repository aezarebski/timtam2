package beast.evolution.speciation;

import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.evolution.tree.birthdeath.BackwardsPointProcess;
import beast.evolution.tree.birthdeath.BackwardsSchedule;
import beast.util.TreeParser;
import org.junit.Test;

import java.util.function.BiPredicate;
import java.util.function.DoubleFunction;

import static org.junit.Assert.assertTrue;

public class TestTimTam {

    private final BiPredicate<Double, Double> approxEqual = (x, y) -> Math.abs(x - y) < 1e-5;
    private final BiPredicate<Double, Double> roughlyEqual = (x, y) -> Math.abs(x - y) < 1e-1;
    // check if within 5% of the second value.
    private final BiPredicate<Double, Double> kindaEqual = (x, y) -> (Math.abs(x - y) / Math.abs(y)) < 5e-2;

    @Test
    public void testLikelihoodCalculationSimple() {
        /**
         * This test draws on a similar one in BDSKY and checks that the TimTam values look similar in a special case.
         *
         * | R0  | lBDSKY | lambda |
         * |-----+--------+--------|
         * | 1.5 | -26.1  |  2.25  |
         * | 1.6 | -27.4  |  2.40  |
         * | 1.7 | -28.8  |  2.55  |
         * | 1.8 | -30.2  |  2.70  |
          */

        TimTam tt =  new TimTam();

        Tree tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
        tt.setInputValue("tree", tree);
        tt.setInputValue("rootLength", new RealParameter("5.0"));
        tt.setInputValue("mu", new RealParameter("1.05"));
        tt.setInputValue("psi", new RealParameter("0.45"));

        tt.setInputValue("lambda", new RealParameter("2.25"));
        tt.initAndValidate();
        assertTrue(kindaEqual.test(-26.1 - 2, tt.calculateLogP()));

        tt.setInputValue("lambda", new RealParameter("2.40"));
        tt.initAndValidate();
        assertTrue(kindaEqual.test(-27.4 - 2, tt.calculateLogP()));

        tt.setInputValue("lambda", new RealParameter("2.55"));
        tt.initAndValidate();
        assertTrue(kindaEqual.test(-28.8 - 2, tt.calculateLogP()));

        tt.setInputValue("lambda", new RealParameter("2.70"));
        tt.initAndValidate();
        assertTrue(kindaEqual.test(-30.2 - 2, tt.calculateLogP()));
    }

    @Test
    public void testLikelihood() {

        RealParameter birthRate = new RealParameter("2.0");
        RealParameter deathRate = new RealParameter("1.0");
        RealParameter samplingRate = new RealParameter("0.3");
        RealParameter rhoProb = new RealParameter("0.5");
        RealParameter occurrenceRate = new RealParameter("0.6");

        RealParameter rootLength = new RealParameter("1.0");
        Tree tree = new TreeParser("(((1:3,2:1):1,3:4):2,4:6);",false);

        BackwardsPointProcess points = new BackwardsPointProcess();
        points.initByName("value", "5.0 1.0");

        BackwardsSchedule catastropheTimes = new BackwardsSchedule();
        catastropheTimes.initByName("value", "0.0");

        TimTam tt = new TimTam(birthRate, deathRate, samplingRate, rhoProb, occurrenceRate, rootLength, catastropheTimes, points);
        tt.setInputValue("tree", tree);

        double fx, fxh, h, fxDash;
        h = 1e-6;
        fxh = tt.p0(1.0, 1.0 - h);
        fx = tt.p0(1.0, 1.0);
        fxDash = tt.lnP0Dash1(1.0);
        assertTrue(approxEqual.test(Math.log((fx - fxh) / h), fxDash));

        assertTrue(
                roughlyEqual.test(
                        -47.0,
                        tt.calculateLogP()
                )
        );
    }

    @Test
    public void testNegativeBinomial() {

        TimTam tt = new TimTam();

        TimTam.NegativeBinomial myNB = tt.getNewNegativeBinomial();
        myNB.setZero();

        assertTrue(
                approxEqual.test(
                        0.0,
                        TimTam.NegativeBinomial.lnPochhammer(1.5, 0)
                )
        );
        assertTrue(
                approxEqual.test(
                        Math.log(1.5),
                        TimTam.NegativeBinomial.lnPochhammer(1.5, 1)
                )
        );
        assertTrue(
                approxEqual.test(
                        Math.log(1.5 * 2.5),
                        TimTam.NegativeBinomial.lnPochhammer(1.5, 2)
                )
        );

        myNB.setLnPAndLnR(Math.log(0.3), Math.log(5));

        assertTrue(approxEqual.test(myNB.getLnMean(),0.7621401));
        assertTrue(approxEqual.test(myNB.getLnVariance(),1.118815));


        DoubleFunction<Double> f = (double z) -> myNB.lnPGF(z, 5, 0.3);
        DoubleFunction<Double> fDash = (double z) -> myNB.lnPGFDash(1, z, 5, 0.3);
        double z = 0.4;
        double h = 1e-6;
        double fFD = Math.log((Math.exp(f.apply(z + h)) - Math.exp(f.apply(z))) / h);

        assertTrue(approxEqual.test(f.apply(z), -1.144208));
        assertTrue(f.apply(z+h) > -1.144208);
        assertTrue(approxEqual.test(fDash.apply(z), fFD));

    }

}
