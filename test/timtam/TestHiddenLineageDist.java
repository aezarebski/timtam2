package timtam;

import org.junit.Test;

import java.util.function.DoubleFunction;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class TestHiddenLineageDist {

    HiddenLineageDist hld = new HiddenLineageDist();

    @Test
    public void testLnPMF() {

        double eps = 1e-7;

        hld.setIsDegenerate(0);

        assertEquals(hld.lnPMF(0), 0.0, eps);
        assertEquals(hld.lnPMF(1), Double.NEGATIVE_INFINITY, eps);
        assertEquals(hld.lnPMF(2), Double.NEGATIVE_INFINITY, eps);
        assertEquals(hld.lnPMF(5), Double.NEGATIVE_INFINITY, eps);
        assertEquals(hld.lnPMF(10), Double.NEGATIVE_INFINITY, eps);

        hld.setIsDegenerate(2);

        assertEquals(hld.lnPMF(0), Double.NEGATIVE_INFINITY, eps);
        assertEquals(hld.lnPMF(1), Double.NEGATIVE_INFINITY, eps);
        assertEquals(hld.lnPMF(2), 0.0, eps);
        assertEquals(hld.lnPMF(5), Double.NEGATIVE_INFINITY, eps);
        assertEquals(hld.lnPMF(10), Double.NEGATIVE_INFINITY, eps);

        hld.setLnPAndLnR(Math.log(0.3), Math.log(5));

        // > options(digits=15)
        // > paste0(as.character(dnbinom(x = c(0,1,2,5,10), size = 5, prob = 1 - 0.3, log = TRUE)), collapse= ",")
        // [1] "-1.78337471969366,-1.3779096115855,-1.48327012724332,-2.96695683437186,-6.9143479836378"
        double[] rTrueVals = new double[] {
            -1.78337471969366,
            -1.3779096115855,
            -1.48327012724332,
            -2.96695683437186,
            -6.9143479836378
        };

        assertEquals(hld.lnPMF(0), rTrueVals[0], eps);
        assertEquals(hld.lnPMF(1), rTrueVals[1], eps);
        assertEquals(hld.lnPMF(2), rTrueVals[2], eps);
        assertEquals(hld.lnPMF(5), rTrueVals[3], eps);
        assertEquals(hld.lnPMF(10), rTrueVals[4], eps);

        hld.setLnMeanAndLnVariance(Math.log(5), Math.log(20));

        assertEquals(hld.getLnP(), Math.log(0.75), eps) ;
        assertEquals(hld.getLnR(), Math.log(1.6666666666), eps) ;

        //        (%i1)	solve([m = n * (1-p) / p, v = n * (1-p)  /(p^2)], [n, p]);
        //        (%o1)	[[n=m^2/(v-m),p=m/v],[n=0,p=0]]

        // > (\(m,v) c(m^2 / (v - m), m / v))(5, 20)
        // [1] 1.66666666666667 0.25000000000000
        // > paste0(as.character(dnbinom(x = c(0,1,2,5,10), size = 1.6666666, prob = 0.25, log = TRUE)), collapse= ",")
        // [1] "-2.31049050944686,-2.08734699813265,-2.08734702313265,-2.47040876647682,-3.49644668161342"
        rTrueVals = new double[] {
            -2.31049059762452,
            -2.08734704931031,
            -2.08734705081031,
            -2.47040875134546,
            -3.49644662590101
        };

        assertEquals(hld.lnPMF(0), rTrueVals[0], eps);
        assertEquals(hld.lnPMF(1), rTrueVals[1], eps);
        assertEquals(hld.lnPMF(2), rTrueVals[2], eps);
        assertEquals(hld.lnPMF(5), rTrueVals[3], eps);
        assertEquals(hld.lnPMF(10), rTrueVals[4], eps);
    }

    @Test
    public void testLnPGF() {
        double[][] testData;
        double eps = 1e-5;

        // PGF is z -> 1
        hld.setIsDegenerate(0);
        testData = new double[][] {
            {1,0.1},{1,0.2},{1,0.5},{1,1.0}
        };
        for (double[] x : testData) {
            assertEquals(Math.log(x[0]), hld.lnPGF(x[1]), eps);
        }

        // PGF is z -> z.
        hld.setIsDegenerate(1);
        testData = new double[][] {
            {0.1,0.1},{0.2,0.2},{0.5,0.5},{1.0,1.0}
        };
        for (double[] x : testData) {
            assertEquals(Math.log(x[0]), hld.lnPGF(x[1]), eps);
        }

        // PGF is z -> z^2.
        hld.setIsDegenerate(2);
        testData = new double[][] {
            {0.01,0.1},{0.04,0.2},{0.25,0.5},{1.0,1.0}
        };
        for (double[] x : testData) {
            assertEquals(Math.log(x[0]), hld.lnPGF(x[1]), eps);
        }

        // PGF is z -> z^100.
        hld.setIsDegenerate(100);
        testData = new double[][] {
            {Math.pow(0.1,100),0.1},
            {Math.pow(0.2,100),0.2},
            {Math.pow(0.5,100),0.5},
            {Math.pow(1.0,100),1.0}
        };
        for (double[] x : testData) {
            assertEquals(Math.log(x[0]), hld.lnPGF(x[1]), eps);
        }

        // PGF is z -> ((1 - p)/(1 - p*z))^r

        // # p = 0.1, r = 4
        // > lnPGF <- \(p,r,z) r * (log(1 - p) - log(1 - p * z))
        // > jstr <- \(xs) xs |> as.character() |> paste0(collapse=",")
        // > foobar <- \(z) sprintf("{%f,%f}", foo(0.1, 4, z), z)
        // > c(0.1,0.2,0.5,1.0) |> sapply(FUN = foobar) |> jstr()
        // [1] "{-0.381241,0.100000},{-0.340631,0.200000},{-0.216269,0.500000},{0.000000,1.000000}"
        testData = new double[][] {
            {-0.381241,0.100000},
            {-0.340631,0.200000},
            {-0.216269,0.500000},
            {0.000000,1.000000}
        };
        hld.setLnPAndLnR(Math.log(0.1), Math.log(4.0));
        for (double[] x : testData) {
            assertEquals(x[0], hld.lnPGF(x[1]), eps);
        }

        // # p = 0.2, r = 100
        // > lnPGF <- \(p,r,z) r * (log(1 - p) - log(1 - p * z))
        // > jstr <- \(xs) xs |> as.character() |> paste0(collapse=",")
        // > foobar <- \(z) sprintf("{%f,%f}", foo(0.2, 100, z), z)
        // > c(0.1,0.2,0.5,1.0) |> sapply(FUN = foobar) |> jstr()
        // [1] "{-20.294084,0.100000},{-18.232156,0.200000},{-11.778304,0.500000},{0.000000,1.000000}"
        testData = new double[][] {
            {-20.294084,0.100000},
            {-18.232156,0.200000},
            {-11.778304,0.500000},
            {0.000000,1.000000}
        };
        hld.setLnPAndLnR(Math.log(0.2), Math.log(100.0));
        for (double[] x : testData) {
            assertEquals(x[0], hld.lnPGF(x[1]), eps);
        }
    }

    @Test
    public void testLnPGFDash() {

        double z, h, fFD;
        DoubleFunction<Double> f, fDash;
        z = 0.4;
        h = 1e-8;

        // Test NB case
        hld.setLnPAndLnR(Math.log(0.3), Math.log(5));
        assertEquals(hld.getLnMean(),0.7621401, 1e-5);
        assertEquals(hld.getLnVariance(),1.118815, 1e-5);
        f = (double x) -> hld.lnPGF(x, 5, 0.3);
        fDash = (double x) -> hld.lnPGFDash(1, x, 5, 0.3);
        fFD = Math.log((Math.exp(f.apply(z + h)) - Math.exp(f.apply(z))) / h);
        assertEquals(f.apply(z), -1.144208, 1e-5);
        assertTrue(f.apply(z+h) > -1.144208);
        assertEquals(fDash.apply(z), fFD, 1e-5);

        // Test zero case
        hld.setIsDegenerate(0);
        assertEquals(hld.getLnMean(), Double.NEGATIVE_INFINITY, 1.0);
        assertEquals(hld.getLnVariance(), Double.NEGATIVE_INFINITY, 1.0);
        f = (double x) -> hld.lnPGF(x);
        fDash = (double x) -> hld.lnPGFDash(1, x);
        fFD = Math.log((Math.exp(f.apply(z + h)) - Math.exp(f.apply(z))) / h);
        assertEquals(f.apply(z), 0.0, 1e-5);
        assertEquals(fDash.apply(z), fFD, 1e-5);

        // Test non-zero degenerate case
        int n = 100;
        hld.setIsDegenerate(n);
        assertEquals(hld.getLnMean(), Math.log(n), 1e-5);
        assertEquals(hld.getLnVariance(), Double.NEGATIVE_INFINITY, 1.0);
        f = (double x) -> hld.lnPGF(x);
        fDash = (double x) -> hld.lnPGFDash(1, x);
        fFD = Math.log((Math.exp(f.apply(z + h)) - Math.exp(f.apply(z))) / h);
        assertEquals(f.apply(z), n * Math.log(z), 1e-5);
        assertEquals(fDash.apply(z), fFD, 1e-5);

    }
}
