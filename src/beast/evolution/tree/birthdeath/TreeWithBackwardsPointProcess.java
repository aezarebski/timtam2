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
public class TreeWithBackwardsPointProcess extends CalculationNode {

    final public Input<RealParameter> rootLengthInput = new Input<>("rootLength", "the time between the origin and the MRCA of the tree", Input.Validate.REQUIRED);
    final public Input<Tree> treeInput = new Input<>("tree", "the tree", Input.Validate.REQUIRED);
    final public Input<BackwardsPointProcess> bwdPointsInput = new Input<>("bwdPoints", "the points in the point process", Input.Validate.REQUIRED);
    final public Input<BackwardsSchedule> bwdCatastropheTimesInput = new Input<>("bwdCatastropheTimes", "the times at which there was a catastrophe", Input.Validate.OPTIONAL);

    private Tree tree;
    private BackwardsSchedule bwdCatastropheTraits;
    private double rootLength;

    public TreeWithBackwardsPointProcess() {
        super();
    }

    public TreeWithBackwardsPointProcess(RealParameter rootLength,
                                         Tree tree,
                                         BackwardsPointProcess bwdPoints,
                                         BackwardsSchedule bwdCatastropheTimes) {
        init(rootLength, tree, bwdPoints, bwdCatastropheTimes);
    }

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        rootLength = rootLengthInput.get().getDoubleValues()[0];
        bwdCatastropheTraits = bwdCatastropheTimesInput.get();

        calculateIntervals();
        intervalsKnown = false;
    }

    /**
     * Stolen from coalescent.TreeIntervals because it is protected there... this function mutates the times array and
     * returns nothing. Note that the times are go backwards from the present. The very last time goes all the way back
     * to the origin.
     *
     * @param tree
     * @param rootLength
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
        // we need to reverse and shift the times in the tree so it uses forward-time starting from the origin. This
        // could probably be improved to avoid the call to the max method.
        double maxTreeNodeTime = Arrays.stream(treeNodeTimes).max().getAsDouble();
        for (int i = 0; i < treeNodeCount; i++) {
            treeNodeTimes[i] = maxTreeNodeTime - treeNodeTimes[i];
        }
        int[] treeJxs = new int[treeNodeCount];
        HeapSort.sort(treeNodeTimes, treeJxs);

        int numCatastrophes, totalInCatastrophes;
        List<Double> fwdCatastropheTimes;
        int[] catastropheSizes;
        if (bwdCatastropheTraits != null) {
            // we make a list of the (forward) times at which there was a catastrophe so that we can collect this
            // information out of the tree later.
            fwdCatastropheTimes = bwdCatastropheTraits.valuesInput.get().stream().mapToDouble(x -> maxTreeNodeTime - x).sorted().boxed().toList();
            catastropheSizes = new int[fwdCatastropheTimes.size()];
            Arrays.fill(catastropheSizes, 0);
            measureCatastrophes(tree, maxTreeNodeTime, fwdCatastropheTimes, catastropheSizes);
            numCatastrophes = catastropheSizes.length;
            totalInCatastrophes = Arrays.stream(catastropheSizes).sum();
        } else {
            // the catastrophe times are null which is used as a signal that there were no catastrophes.
            fwdCatastropheTimes = new ArrayList<>();
            numCatastrophes = 0;
            totalInCatastrophes = 0;
            catastropheSizes = new int[0];
        }

        double[] fwdPointTimes = bwdPointsInput.get().valuesInput.get().stream().mapToDouble(x -> maxTreeNodeTime - x).sorted().toArray();
        int pointCount = fwdPointTimes.length;

        // The number of intervals needs to account for catastrophes where there are no sequences collected and
        // catastrophes where multiple leaves correspond to a single interval.
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
            pointET = (pointJx < fwdPointTimes.length) ? fwdPointTimes[pointJx] : Double.POSITIVE_INFINITY;
            // we keep track of the time of the next catastrophe so that it is easier to check how many sequences where sampled and when it occurred.
            catastET = (catastJx < fwdCatastropheTimes.size()) ? fwdCatastropheTimes.get(catastJx) : Double.POSITIVE_INFINITY;

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
     * @param maxTime
     * @param fwdCatastropheTimes
     * @param catastropheSizes
     *
     * @see TreeWithBackwardsPointProcess#collectTimes
     */
    private void measureCatastrophes(Tree tree, double maxTime, List<Double> fwdCatastropheTimes, int[] catastropheSizes) {
        Node[] nodes = tree.getNodesAsArray();
        double nodeTime;

        for (Node node : nodes) {
            if (node.isLeaf()) {
                nodeTime = maxTime - node.getHeight();
                for (int j = 0; j < fwdCatastropheTimes.size(); j++) {
                    if (Math.abs(fwdCatastropheTimes.get(j) - nodeTime) < 1e-8) {
                        catastropheSizes[j] += 1;
                    }
                }
            }
        }

        if (Arrays.stream(catastropheSizes).sum() == 0) {
            StringBuilder nodeTimesSB = new StringBuilder();
            for (Node node : nodes) {
                if (node.isLeaf()) {
                    nodeTimesSB.append("\n\t\t").append(maxTime - node.getHeight());
                }
            }
            throw new RuntimeException(
                    "\nIt appears that no samples were obtained in any of the catastrophies, please " +
                            "\ndouble check that the timing of these has been specified correctly in the XML. " +
                            "\nFor debugging purposes here are some variables:" +
                            "\n\tThe current catastrophe times are " + fwdCatastropheTimes +
                            "\n\tThe maxTime is " + maxTime +
                            "\n\tThe node times are " + nodeTimesSB);
        }
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