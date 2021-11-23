package beast.evolution.tree.birthdeath;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.RandomTree;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.ConstantPopulation;
import beast.math.distributions.MRCAPrior;
import beast.util.TreeParser;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;
import java.util.function.BiPredicate;

public class TestTreeWithBackwardsPointProcess {

    private final BiPredicate<Double, Double> approxEqual = (x, y) -> Math.abs(x - y) < 1e-5;

    @Test
    public void estimatedTreeEndingOnOccurrence() {

        ConstantPopulation myPopModel = new ConstantPopulation();
        myPopModel.initByName("popSize", "1.0");

        List<Sequence> mySeqs = new ArrayList<>();
        mySeqs.add(new Sequence("t0_8.0", "?"));
        mySeqs.add(new Sequence("t1_9.0", "?"));
        Alignment myTaxa = new Alignment(mySeqs, "nucleotide");
        TaxonSet myTaxonSet = new TaxonSet(myTaxa);

        TraitSet myTraitSet = new TraitSet();
        myTraitSet.initByName(
                "traitname", "date",
                "value", "t0_8.0=8.0,t1_9.0=9.0",
                "taxa", myTaxonSet
        );

        RandomTree myTree;
        myTree = new RandomTree();
        myTree.initByName(
                "taxonset", myTaxonSet,
                "trait", myTraitSet,
                "populationModel", myPopModel);

        double rootLengthDouble = 1.0;
        RealParameter rootLength = new RealParameter("1.0");

        // This loop is included to ensure that we get an origin that occurs before the first occurrence event.
        double mRCAHeight;
        mRCAHeight = myTree.getRoot().getHeight();
        while (mRCAHeight + rootLengthDouble < 4.0) {
            myTree = new RandomTree();
            myTree.initByName(
                    "taxonset", myTaxonSet,
                    "trait", myTraitSet,
                    "populationModel", myPopModel);
            mRCAHeight = myTree.getRoot().getHeight();
        }
        assertTrue(mRCAHeight + rootLengthDouble > 4.0);

        BackwardsPointProcess points = new BackwardsPointProcess();
        points.initByName("value", "4.0 -0.5");

        TreeWithBackwardsPointProcess tpp = new TreeWithBackwardsPointProcess(
                rootLength,
                myTree,
                points,
                null);

        // There are two occurrences, two leaves and one internal nodes. The "present" is the time of the last tip.
        assertEquals(5, tpp.getIntervalCount());

        double hiddenDelay;
        if (Objects.equals(tpp.getIntervalType(0).toString(), "birth")) {
            assertEquals(tpp.getIntervalType(0).toString(), "birth");
            assertTrue(approxEqual.test(tpp.getIntervalDuration(0), 1.0));
            assertEquals(tpp.getIntervalType(1).toString(), "occurrence");
            assertTrue(tpp.getIntervalDuration(1) > 0.0);
            hiddenDelay = tpp.getIntervalDuration(1) + 1.0;
        } else {
            assertEquals(tpp.getIntervalType(0).toString(), "occurrence");
            assertTrue(tpp.getIntervalDuration(0) > 0);
            assertEquals(tpp.getIntervalType(1).toString(), "birth");
            assertTrue(tpp.getIntervalDuration(1) < 1.0);
            hiddenDelay = tpp.getIntervalDuration(0) + tpp.getIntervalDuration(1);
        }

        assertEquals(tpp.getIntervalType(2).toString(), "sample");
        assertTrue(approxEqual.test(
                hiddenDelay + tpp.getIntervalDuration(2),
                mRCAHeight + rootLengthDouble - 1.0));


        assertEquals(tpp.getIntervalType(3).toString(), "sample");
        assertTrue(approxEqual.test(tpp.getIntervalDuration(3), 1.0));

        assertEquals(tpp.getIntervalType(4).toString(), "occurrence");
        assertTrue(approxEqual.test(tpp.getIntervalDuration(4), 0.5));

    }

    @Test
    public void estimatedTreeEndingOnSample() {

        ConstantPopulation myPopModel = new ConstantPopulation();
        myPopModel.initByName("popSize", "1.0");

        List<Sequence> mySeqs = new ArrayList<>();
        mySeqs.add(new Sequence("t0_8.0", "?"));
        mySeqs.add(new Sequence("t1_9.5", "?"));
        Alignment myTaxa = new Alignment(mySeqs, "nucleotide");
        TaxonSet myTaxonSet = new TaxonSet(myTaxa);

        TraitSet myTraitSet = new TraitSet();
        myTraitSet.initByName(
                "traitname", "date",
                "value", "t0_8.0=8.0,t1_9.5=9.5",
                "taxa", myTaxonSet
        );

        RandomTree myTree;
        myTree = new RandomTree();
        myTree.initByName(
                "taxonset", myTaxonSet,
                "trait", myTraitSet,
                "populationModel", myPopModel);

        double rootLengthDouble = 1.0;
        RealParameter rootLength = new RealParameter("1.0");

        // This loop is included to ensure that we get an origin that occurs before the first occurrence event.
        double mRCAHeight;
        mRCAHeight = myTree.getRoot().getHeight();
        while (mRCAHeight + rootLengthDouble < 4.5) {
            myTree = new RandomTree();
            myTree.initByName(
                    "taxonset", myTaxonSet,
                    "trait", myTraitSet,
                    "populationModel", myPopModel);
            mRCAHeight = myTree.getRoot().getHeight();
        }
        assertTrue(mRCAHeight + rootLengthDouble > 4.5);


        BackwardsPointProcess points = new BackwardsPointProcess();
        points.initByName("value", "4.5 0.5");

        TreeWithBackwardsPointProcess tpp = new TreeWithBackwardsPointProcess(
                rootLength,
                myTree,
                points,
                null);

        // There are two occurrences, two leaves and one internal nodes. The "present" is the time of the last tip.
        assertEquals(5, tpp.getIntervalCount());

        double hiddenDelay;
        if (Objects.equals(tpp.getIntervalType(0).toString(), "birth")) {
            assertEquals(tpp.getIntervalType(0).toString(), "birth");
            assertTrue(approxEqual.test(tpp.getIntervalDuration(0), 1.0));
            assertEquals(tpp.getIntervalType(1).toString(), "occurrence");
            assertTrue(tpp.getIntervalDuration(1) > 0.0);
            hiddenDelay = tpp.getIntervalDuration(1) + 1.0;
        } else {
            assertEquals(tpp.getIntervalType(0).toString(), "occurrence");
            assertTrue(tpp.getIntervalDuration(0) > 0);
            assertEquals(tpp.getIntervalType(1).toString(), "birth");
            assertTrue(tpp.getIntervalDuration(1) < 1.0);
            hiddenDelay = tpp.getIntervalDuration(0) + tpp.getIntervalDuration(1);
        }

        assertEquals(tpp.getIntervalType(2).toString(), "sample");
        assertTrue(approxEqual.test(
                hiddenDelay + tpp.getIntervalDuration(2),
                mRCAHeight + rootLengthDouble - 1.5));


        assertEquals(tpp.getIntervalType(3).toString(), "occurrence");
        assertTrue(approxEqual.test(tpp.getIntervalDuration(3), 1.0));

        assertEquals(tpp.getIntervalType(4).toString(), "sample");
        assertTrue(approxEqual.test(tpp.getIntervalDuration(4), 0.5));

    }

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
