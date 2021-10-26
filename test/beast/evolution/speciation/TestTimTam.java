package beast.evolution.speciation;

import org.junit.Test;

import java.util.function.BiPredicate;
import java.util.function.DoubleFunction;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class TestTimTam {

    private BiPredicate<Double, Double> approxEqual = (x, y) -> Math.abs(x - y) < 1e-5;

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
