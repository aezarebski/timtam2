import beast.core.parameter.RealParameter;
import beast.evolution.speciation.TimTam;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;

public class Demo {
    public static void main(String[] args) {

        TimTam tt = new TimTam();

        // the false argument is so that the tips are not adjusted to make the tree ultrametric.
        Tree tree = new TreeParser("((1: 4.5, 2: 4.5):1,3:5.5);",false);

        tt.setInputValue("tree", tree);
        tt.setInputValue("lambda", new RealParameter("1.5"));
        tt.setInputValue("mu", new RealParameter("0.3"));
        tt.setInputValue("psi", new RealParameter("0.3"));
        tt.setInputValue("p", new RealParameter("0.3"));
        tt.setInputValue("hasFinalSample", false);
        tt.setInputValue("origin", new RealParameter("6.5"));

        tt.initAndValidate();
        double foo = tt.calculateLogP();
        System.out.println("Hello World" + foo);
    }
}
