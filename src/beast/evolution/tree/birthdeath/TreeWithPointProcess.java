package beast.evolution.tree.birthdeath;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.HeapSort;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.OptionalInt;

@Description("Extracts intervals from a tree when there is an additional point-process associated with it.")
public class TreeWithPointProcess extends CalculationNode {

    final public Input<RealParameter> rootLengthInput = new Input<>("rootLength", "the time between the origin and the MRCA of the tree", Input.Validate.REQUIRED);
    final public Input<Tree> treeInput = new Input<>("tree", "the tree", Input.Validate.REQUIRED);
    final public Input<PointProcess> pointsInput = new Input<>("points", "the points in the point process", Input.Validate.REQUIRED);
    final public Input<Schedule> catastropheTimesInput = new Input<>("catastropheTimes", "the times at which there was a catastrophe", Input.Validate.OPTIONAL);

    private Tree tree;
    private Schedule catastropheTraits;
    private double rootLength;
    private double[] pointTimes;
    private int pointCount;

    public TreeWithPointProcess() {
        super();
    }

    public TreeWithPointProcess(RealParameter rootLength, Tree tree, PointProcess points, Schedule catastropheTimes) {
        init(rootLength, tree, points, catastropheTimes);
    }

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        catastropheTraits = catastropheTimesInput.get();
        rootLength = rootLengthInput.get().getDoubleValues()[0];

        PointProcess pointTraits = pointsInput.get();
        pointTimes = pointTraits.valuesInput.get().stream().sorted().mapToDouble(Double::doubleValue).toArray();
        pointCount = pointTimes.length;
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

        int numCatastrophes, totalInCatastrophes;
        List<Double> catastropheTimes;
        int[] catastropheSizes;
        if (catastropheTraits != null) {
            // we make a list of the (forward) times at which there was a catastrophe so that we can collect this information out of the tree later.
            catastropheTimes = catastropheTraits.valuesInput.get().stream().mapToDouble(Double::doubleValue).sorted().boxed().toList();
            catastropheSizes = new int[catastropheTimes.size()];
            Arrays.fill(catastropheSizes, 0);
            measureCatastrophes(tree, maxTreeNodeTime, catastropheTimes, catastropheSizes);
            numCatastrophes = catastropheSizes.length;
            totalInCatastrophes = Arrays.stream(catastropheSizes).sum();
        } else {
            // the catastrophe times are null which us used as a signal that there were no catastrophes.
            catastropheTimes = new ArrayList<>();
            numCatastrophes = 0;
            totalInCatastrophes = 0;
            catastropheSizes = new int[0];
        }


        // The number of intervals needs to account for catastrophes where there are no sequences collected and catastrophes where multiple leaves correspond to a single interval.
        intervalCount = pointCount + treeNodeCount - totalInCatastrophes + numCatastrophes;
        intervals = new double[intervalCount];
        intervalTypes = new EventType[intervalCount];

        int treeJx = 0;
        int pointJx = 0;
        int catastJx = 0;
        double treeET, pointET, catastET;
        double currTime = 0.0; // TODO we should not really be hardcoding that the origin is zero....
        int intIx = 0;
        while (intIx < intervalCount) {
            // if the tree events have been exhausted set the next event time to positive
            // infinity so that it cannot appear as the next observation. Do the same the point process.
            treeET = (treeJx < treeJxs.length) ? treeNodeTimes[treeJxs[treeJx]] : Double.POSITIVE_INFINITY;
            pointET = (pointJx < pointTimes.length) ? pointTimes[pointJx] : Double.POSITIVE_INFINITY;
            // we keep track of the time of the next catastrophe so that it is easier to check how many sequences where sampled and when it occurred.
            catastET = (catastJx < catastropheTimes.size()) ? catastropheTimes.get(catastJx) : Double.POSITIVE_INFINITY;

            if (treeET < pointET) {
                intervals[intIx] = treeET - currTime;
                currTime = treeET;
                if (Math.abs(treeET - catastET) > 1e-7) {
                    if (treeNodeOutdegree[treeJxs[treeJx]] == 2) {
                        intervalTypes[intIx] = new EventType("birth", OptionalInt.empty());
                    } else if (treeNodeOutdegree[treeJxs[treeJx]] == 0) {
                        intervalTypes[intIx] = new EventType("sample", OptionalInt.empty());
                    } else {
                        throw new IllegalArgumentException("Non-binary tree provided to TreeWithPointProcess");
                    }
                    treeJx++;
                } else {
                    intervalTypes[intIx] = new EventType("catastrophe", OptionalInt.of(catastropheSizes[catastJx]));
                    catastJx++;
                    while (Math.abs(treeET - catastET) < 1e-7) {
                        treeJx++;
                        treeET = (treeJx < treeJxs.length) ? treeNodeTimes[treeJxs[treeJx]] : Double.POSITIVE_INFINITY;
                    }
                }
            } else if (treeET > pointET) {
                intervals[intIx] = pointET - currTime;
                currTime = pointET;
                intervalTypes[intIx] = new EventType("occurrence", OptionalInt.empty());
                pointJx++;
            } else {
                throw new IllegalArgumentException("It appears that a tree event and a point process event occurred at the same time:\ntree time " + treeET + ", point time " + pointET);
            }
            intIx++;
        }

        intervalsKnown = true;
    }

    /**
     * Look at the tree and detect catastrophes and mutate the array of sizes to contain these counts.
     *
     * @param tree
     * @param catastropheTimes
     * @param catastropheSizes
     *
     * @see TreeWithPointProcess#collectTimes
     */
    private void measureCatastrophes(Tree tree, double maxTime, List<Double> catastropheTimes, int[] catastropheSizes) {
        Node[] nodes = tree.getNodesAsArray();
        double nodeTime;

        for (Node node : nodes) {
            if (node.isLeaf()) {
                nodeTime = maxTime - node.getHeight();
                for (int j = 0; j < catastropheTimes.size(); j++) {
                    if (Math.abs(catastropheTimes.get(j) - nodeTime) < 1e-8) {
                        catastropheSizes[j] += 1;
                    }
                }
            }
        }

        if (Arrays.stream(catastropheSizes).sum() == 0) throw new RuntimeException(
                "\nIt appears that no samples were obtained in any of the catastrophies, please " +
                "\ndouble check that the timing of these has been specified correctly in the XML. " +
                "\nThe current catastrophe times are " + catastropheTimes);
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
        if (i < 0 || i >= intervalCount) throw new IllegalArgumentException("intervalCount is " + intervalCount + " but interval " + i + " requested.");
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