package beast.evolution.tree.birthdeath;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.util.HeapSort;

import java.util.Arrays;

@Description("Extracts intervals from a tree when there is an additional point-process associated with it.")
public class TreeWithPointProcess extends CalculationNode {

    final public Input<RealParameter> originInput = new Input<>("origin", "the origin time", Input.Validate.REQUIRED);
    final public Input<Tree> treeInput = new Input<>("tree", "the tree", Input.Validate.REQUIRED);
    final public Input<TraitSet> pointsInput = new Input<>("points", "the points in the point process", Input.Validate.REQUIRED);

    public TreeWithPointProcess() {
        super();
    }

    public TreeWithPointProcess(RealParameter origin, Tree tree, TraitSet points) {
        init(origin, tree, points);
    }

    @Override
    public void initAndValidate() {
        calculateIntervals();
        intervalsKnown = false;
    }

    /**
     * Stolen from coalescent.TreeIntervals because it is protected there...
     *
     * @param tree
     * @param times
     * @param childCounts
     */
    protected static void collectTimes(Tree tree, double[] times, int[] childCounts) {
        Node[] nodes = tree.getNodesAsArray();
        for (int i = 0; i < nodes.length; i++) {
            Node node = nodes[i];
            times[i] = node.getHeight();
            childCounts[i] = node.isLeaf() ? 0 : 2;
        }
    }

    /**
     * Calculate the intervals.
     */
    protected void calculateIntervals() {
        Tree tree = treeInput.get();
        TraitSet pointTraits = pointsInput.get();

        final int treeNodeCount = tree.getNodeCount();
        double[] treeNodeTimes = new double[treeNodeCount];
        int[] treeNodeOutdegree = new int[treeNodeCount];
        collectTimes(tree, treeNodeTimes, treeNodeOutdegree);
        // we need to reverse and shift the times in the tree so it uses forward-time
        // starting from the origin.
        // TODO ensure that this is handling the difference between TMRCA and the origin time correctly.
        double maxTreeNodeTime = Arrays.stream(treeNodeTimes).max().getAsDouble();
        for (int i = 0; i < treeNodeCount; i++) {
            treeNodeTimes[i] = maxTreeNodeTime - treeNodeTimes[i];
        }
        int[] treeJxs = new int[treeNodeCount];
        HeapSort.sort(treeNodeTimes, treeJxs);

        final int pointCount = pointTraits.taxaInput.get().getTaxonCount();
        final double[] pointTimes = pointTraits.taxaInput.get().asStringList().stream().mapToDouble((tn) -> pointTraits.getValue(tn)).sorted().toArray();

        intervalCount = treeNodeCount + pointCount;
        intervals = new double[intervalCount];
        intervalTypes = new EventType[intervalCount];

        int treeJx = 0;
        double treeET;
        int pointJx = 0;
        double pointET;
        double currTime = originInput.get().getValue();
        int intIx = 0;
        while (intIx < intervalCount) {
            // if the tree events have been exhausted set the next event time to positive
            // infinity so that it cannot appear as the next observation.
            treeET = (treeJx < treeJxs.length) ? treeNodeTimes[treeJxs[treeJx]] : Double.POSITIVE_INFINITY;
            pointET = (pointJx < pointTimes.length) ? pointTimes[pointJx] : Double.POSITIVE_INFINITY;

            if (treeET < pointET) {
                intervals[intIx] = treeET - currTime;
                currTime = treeET;
                if (treeNodeOutdegree[treeJxs[treeJx]] == 2) {
                    intervalTypes[intIx] = EventType.BIRTH;
                } else if (treeNodeOutdegree[treeJxs[treeJx]] == 0) {
                    intervalTypes[intIx] = EventType.SAMPLE;
                } else {
                    throw new IllegalArgumentException("Non-binary tree provided to TreeWithPointProcess");
                }
                treeJx++;
            } else if (treeET > pointET) {
                intervals[intIx] = pointET - currTime;
                currTime = pointET;
                intervalTypes[intIx] = EventType.OCCURRENCE;
                pointJx++;
            } else {
                throw new IllegalArgumentException("Invalid events given to TreeWithPointProcess");
            }
            intIx++;
        }

        intervalsKnown = true;
    }

    protected int intervalCount;

    protected double[] intervals;

    protected EventType[] intervalTypes;

    protected boolean intervalsKnown = false;


    public int getIntervalCount() {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        return intervalCount;
    }

    public EventType getIntervalType(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        if (i < 0 || i >= intervalCount) throw new IllegalArgumentException();
        return intervalTypes[i];
    }

    /**
     * The time since the last observed event. We adopt the convention that the origin
     * is event number -1 so that events after that start at 0.
     *
     * @param i the event number
     * @return the time between events i-1 and i
     */
    public double getIntervalDuration(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        if (i < 0 || i >= intervalCount) throw new IllegalArgumentException();
        return intervals[i];
    }
}