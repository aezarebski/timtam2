package beast.evolution.tree.birthdeath;

import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.function.BiPredicate;

public class TestTreeWithBackwardsPointProcess {

    private final BiPredicate<Double, Double> approxEqual = (x, y) -> Math.abs(x - y) < 1e-5;

    @Test
    public void canRecoverCatastrophe() {

        RealParameter rootLength = new RealParameter("1.0");
        Tree tree = new TreeParser("(((1:3,2:1):1,3:4):2,4:6);",false);

        BackwardsPointProcess points = new BackwardsPointProcess();
        points.initByName("value", "5.0 1.0");

        BackwardsSchedule catastropheTimes = new BackwardsSchedule();
        catastropheTimes.initByName("value", "0.0");

        TreeWithBackwardsPointProcess tpp = new TreeWithBackwardsPointProcess(rootLength, tree, points, catastropheTimes);

        // There are two occurrences, three leaves and two internal nodes.
        assertEquals(7, tpp.getIntervalCount());

        // the first event is a birth at time 1.0
        assertEquals(tpp.getIntervalType(0).toString(), "birth");
        assertTrue(approxEqual.test(tpp.getIntervalDuration(0), 1.0));

        assertEquals(tpp.getIntervalType(1).toString(), "occurrence");
        assertTrue(approxEqual.test(tpp.getIntervalDuration(1), 1.0));

        assertEquals(tpp.getIntervalType(2).toString(), "birth");
        assertTrue(approxEqual.test(tpp.getIntervalDuration(2), 1.0));

        assertEquals(tpp.getIntervalType(3).toString(), "birth");
        assertTrue(approxEqual.test(tpp.getIntervalDuration(3), 1.0));

        assertEquals(tpp.getIntervalType(4).toString(), "sample");
        assertTrue(approxEqual.test(tpp.getIntervalDuration(4), 1.0));

        assertEquals(tpp.getIntervalType(5).toString(), "occurrence");
        assertTrue(approxEqual.test(tpp.getIntervalDuration(5), 1.0));

        assertEquals(tpp.getIntervalType(6).toString(), "catastrophe");
        assertEquals(tpp.getIntervalType(6).getCount().getAsInt(), 3);
        assertTrue(approxEqual.test(tpp.getIntervalDuration(6), 1.0));
    }

    @Test
    public void canRecoverIntervals() {

        RealParameter rootLength = new RealParameter("1.0");
        Tree tree = new TreeParser("((2:2, 1:1):1,3:4);",false);

        BackwardsPointProcess points = new BackwardsPointProcess();
        points.initByName("value", "3.5 2.5");

        TreeWithBackwardsPointProcess tpp = new TreeWithBackwardsPointProcess(rootLength, tree, points, null);

        // There are two occurrences, three leaves and two internal nodes.
        assertEquals(tpp.getIntervalCount(),2 + 3 + (3 - 1));

        // the first event is a birth at time 1.0
        assertEquals(tpp.getIntervalType(0).toString(), "birth");
        assertTrue(approxEqual.test(tpp.getIntervalDuration(0), 1.0));

        assertEquals(tpp.getIntervalType(1).toString(), "occurrence");
        assertTrue(approxEqual.test(tpp.getIntervalDuration(1), 0.5));

        assertEquals(tpp.getIntervalType(2).toString(), "birth");
        assertTrue(approxEqual.test(tpp.getIntervalDuration(2), 0.5));

        assertEquals(tpp.getIntervalType(3).toString(), "occurrence");
         assertTrue(approxEqual.test(tpp.getIntervalDuration(3), 0.5));

        assertEquals(tpp.getIntervalType(4).toString(), "sample");
        assertTrue(approxEqual.test(tpp.getIntervalDuration(4), 0.5));

        assertEquals(tpp.getIntervalType(5).toString(), "sample");
        assertTrue(approxEqual.test(tpp.getIntervalDuration(5), 1.0));

        assertEquals(tpp.getIntervalType(6).toString(), "sample");
        assertTrue(approxEqual.test(tpp.getIntervalDuration(6), 1.0));
    }
}
