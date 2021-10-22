import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.speciation.TimTam;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.evolution.tree.birthdeath.TreeWithPointProcess;
import beast.evolution.tree.coalescent.TreeIntervals;
import beast.util.TreeParser;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class Demo {

    /** A dummy taxonset.
     *
     * @param Nleaves the number of taxa
     * @return a TaxonSet containing taxa with IDs: t0, t1, ...
     */
    public static TaxonSet dummyTaxonSet(int Nleaves) {
        List<Sequence> seqList = new ArrayList<>();

        for (int i=0; i<Nleaves; i++) {
            String taxonID = "t" + i;
            seqList.add(new Sequence(taxonID, "?"));
        }

        Alignment alignment = new Alignment(seqList, "nucleotide");
        TaxonSet taxonSet = new TaxonSet(alignment);
        return taxonSet;
    }

    public static void main(String[] args) {


        // the false argument is so that the tips are not adjusted to make the tree ultrametric.
        Tree tree = new TreeParser("((1: 4.5, 2: 4.5):1,3:5.5);",false);

//        TimTam tt = new TimTam();
//        tt.setInputValue("tree", tree);
//        tt.setInputValue("lambda", new RealParameter("1.5"));
//        tt.setInputValue("mu", new RealParameter("0.3"));
//        tt.setInputValue("psi", new RealParameter("0.3"));
//        tt.setInputValue("p", new RealParameter("0.3"));
//        tt.setInputValue("hasFinalSample", false);
//        tt.setInputValue("origin", new RealParameter("6.5"));
//        tt.initAndValidate();
//        double foo = tt.calculateLogP();
//        System.out.println("Hello World" + foo);


        RealParameter origin = new RealParameter("0.0");

        // We construct a TraitSet to represent the dates of each of the unsequenced
        // samples because this allows us to create tips without including them as part of
        // the tree but allows for more structure than just a list of numbers.
        TraitSet points = new TraitSet();
        points.initByName("traitname", "point-date",
                "taxa", dummyTaxonSet(2),
                "value", "t0=0.1, t1=4");

        TreeWithPointProcess tpp = new TreeWithPointProcess(origin, tree, points);
        int num_intervals_2 = tpp.getIntervalCount();
        for (int i = 0; i < num_intervals_2; i++) {
            System.out.println("----------------");
            System.out.println(tpp.getIntervalType(i).toString());
            System.out.println(tpp.getIntervalDuration(i));
        }

//        TreeIntervals ti = new TreeIntervals(tree);
//        int num_intervals = ti.getIntervalCount();
//        for (int i = num_intervals-1; i >= 0; i--) {
//            System.out.println("----------------");
//            System.out.println(ti.getIntervalType(i));
//            System.out.println(ti.getInterval(i));
//        }
    }
}
