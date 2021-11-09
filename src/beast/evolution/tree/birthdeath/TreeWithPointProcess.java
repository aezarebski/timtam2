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
import java.util.OptionalInt;

@Description("Extracts intervals from a tree when there is an additional point-process associated with it.")
public class TreeWithPointProcess extends CalculationNode {

    final public Input<RealParameter> rootLengthInput = new Input<>("rootLength", "the time between the origin and the MRCA of the tree", Input.Validate.REQUIRED);
    final public Input<Tree> treeInput = new Input<>("tree", "the tree", Input.Validate.REQUIRED);
    final public Input<TraitSet> pointsInput = new Input<>("points", "the points in the point process", Input.Validate.REQUIRED);

    public TreeWithPointProcess() {
        super();
    }

    public TreeWithPointProcess(RealParameter rootLength, Tree tree, TraitSet points) {
        init(rootLength, tree, points);
    }

    @Override
    public void initAndValidate() {
        calculateIntervals();
        intervalsKnown = false;
    }

    /**
     * Stolen from coalescent.TreeIntervals because it is protected there... this function mutates the times array and
     * returns nothing. Note that the times are go backwards from the present.
     *
     * @param tree
     * @param times
     * @param childCounts
     */
    protected static void collectTimes(Tree tree, double rootLength, double[] times, int[] childCounts) {
        Node[] nodes = tree.getNodesAsArray();
        double maxHeight = 0;
        double currHeight;
        for (int i = 0; i < nodes.length; i++) {
            Node node = nodes[i];
            currHeight = node.getHeight();
            if (currHeight > maxHeight)
                maxHeight = currHeight;
            times[i] = currHeight;
            childCounts[i] = node.isLeaf() ? 0 : 2;
        }
        times[nodes.length] = maxHeight + rootLength;
    }

    /**
     * Calculate the intervals.
     */
    protected void calculateIntervals() {
        Tree tree = treeInput.get();
        TraitSet pointTraits = pointsInput.get();

        double rootLength = rootLengthInput.get().getDoubleValues()[0];

        final int treeNodeCount = tree.getNodeCount();
        double[] treeNodeTimes = new double[treeNodeCount+1]; // one extra time for the origin to tMCRA.
        int[] treeNodeOutdegree = new int[treeNodeCount];
        collectTimes(tree, rootLength, treeNodeTimes, treeNodeOutdegree);
        // we need to reverse and shift the times in the tree so it uses forward-time
        // starting from the origin.
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
        double currTime = 0.0; // TODO we should not really be hardcoding that the origin is zero....
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
                    intervalTypes[intIx] = new EventType("birth", OptionalInt.empty());
                } else if (treeNodeOutdegree[treeJxs[treeJx]] == 0) {
                    intervalTypes[intIx] = new EventType("sample", OptionalInt.empty());
                } else {
                    throw new IllegalArgumentException("Non-binary tree provided to TreeWithPointProcess");
                }
                treeJx++;
            } else if (treeET > pointET) {
                intervals[intIx] = pointET - currTime;
                currTime = pointET;
                intervalTypes[intIx] = new EventType("occurrence", OptionalInt.empty());
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

    public double getInterval(int i) {
        if (!intervalsKnown) {
            calculateIntervals();
        }
        if (i < 0 || i >= intervalCount) throw new IllegalArgumentException();
        return intervals[i];
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