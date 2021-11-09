package beast.evolution.speciation;

import beast.core.parameter.RealParameter;
import beast.evolution.branchratemodel.RateStatistic;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import org.junit.Test;

import java.util.function.BiPredicate;
import java.util.function.DoubleFunction;

import static beast.evolution.tree.birthdeath.TestTreeWithPointProcess.dummyTaxonSet;
import static org.junit.Assert.assertTrue;

public class TestTimTam {

    private final BiPredicate<Double, Double> approxEqual = (x, y) -> Math.abs(x - y) < 1e-5;
    private final BiPredicate<Double, Double> roughlyEqual = (x, y) -> Math.abs(x - y) < 1e-1;

    @Test
    public void testLikelihood() {

        RealParameter birthRate = new RealParameter("2.0");
        RealParameter deathRate = new RealParameter("1.0");
        RealParameter samplingRate = new RealParameter("0.3");
        RealParameter rhoProb = new RealParameter("0.5");
        RealParameter occurrenceRate = new RealParameter("0.6");

        RealParameter rootLength = new RealParameter("1.0");
        Tree tree = new TreeParser("(((1:3,2:1):1,3:4):2,4:6);",false);

        // TODO this can be removed if you use this trick: https://github.com/laduplessis/skylinetools/blob/7993f64ca5699ad0da053af2f2a6d8e8a4c45fa8/src/skylinetools/parameter/DateParser.java#L33
        TraitSet points = new TraitSet();
        points.initByName("traitname", "point-date",
                "taxa", dummyTaxonSet(2),
                "value", "t0=2.0, t1=6.0");

        TraitSet catastropheTimes = new TraitSet();
        catastropheTimes.initByName("traitname", "catastrophe-date",
                "taxa", dummyTaxonSet(1),
                "value", "t0=7.0");

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

        TimTam.NegativeBinomial myNB = tt.getNegativeBinomial();


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
