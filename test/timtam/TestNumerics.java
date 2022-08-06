package timtam;

import org.junit.Test;

import static org.junit.Assert.assertEquals;


public class TestNumerics {

    @Test
    public void testLnPochhammer() {
        double[] xs = new double[]{1.5, 2.7};
        double eps = 1e-14;

        for (double x : xs) {
            assertEquals(0.0, Numerics.lnPochhammer(x, 0), eps);
            assertEquals(Math.log(x), Numerics.lnPochhammer(x, 1), eps);
            assertEquals(
                Math.log(x * (x + 1.0)),
                Numerics.lnPochhammer(x, 2),
                eps);
            assertEquals(
                Math.log(x * (x + 1.0) * (x + 2.0)),
                Numerics.lnPochhammer(x, 3),
                eps);
        }
    }

    @Test
    public void testLnChoose() {

        double eps = 1e-10;

        // > lchoose(n = c(5,5.5,6,6.5), k = c(2,3,4,0))
        // [1] 2.30258509299405 2.66982898828201 2.70805020110221 0.00000000000000

        double[][] rResults = new double[][] {
            {5,2,2.30258509299405},
            {5.5,3,2.66982898828201},
            {6,4,2.70805020110221},
            {6.5,0,0.00000000000000}
        };
        for (double[] vals : rResults) {
            assertEquals(
                Numerics.lnChoose(vals[0], (int) vals[1]),
                vals[2],
                eps
            );
        }
    }

    @Test
    public void testLogSumExp() {

        double eps = 1e-10;

        // > log(sum(exp(c(1.17,0.95,0.10,1.58))))
        // [1] 2.46536945456999
        // > log(sum(exp(c(1.17,0.95,0.10))))
        // [1] 1.93338535684413

        double[] xs = {1.17,0.95,0.10,1.58};

        assertEquals(2.46536945456999, Numerics.logSumExp(xs, 4), eps);
        assertEquals(1.93338535684413, Numerics.logSumExp(xs, 3), eps);


        // > log(sum(exp(c(-1.0,0.95,0.10,1.58))))
        // [1] 2.18759119492718
        xs[0] = -1.0;
        assertEquals(2.18759119492718, Numerics.logSumExp(xs), eps);
        // > log(sum(exp(c(1.0,0.95,0.10,1.58))))
        // [1] 2.42162229448561
        xs[0] = 1.0;
        assertEquals(2.42162229448561, Numerics.logSumExp(xs), eps);
        // > log(sum(exp(c(10.0,0.95,0.10,1.58))))
        // [1] 10.0003879051269
        xs[0] = 10.0;
        assertEquals(10.0003879051269, Numerics.logSumExp(xs), eps);
    }
}
