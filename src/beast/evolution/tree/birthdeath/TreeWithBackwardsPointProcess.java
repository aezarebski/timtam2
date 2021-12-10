package beast.evolution.tree.birthdeath;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.HeapSort;

import java.util.Arrays;
import java.util.OptionalInt;

@Description("Extracts intervals from a tree when there is an additional point-process associated with it.")
public class TreeWithBackwardsPointProcess extends CalculationNode {

    final public Input<RealParameter> rootLengthInput = new Input<>("rootLength", "the time between the origin and the MRCA of the tree", Input.Validate.REQUIRED);
    final public Input<Tree> treeInput = new Input<>("tree", "the tree", Input.Validate.REQUIRED);
    final public Input<BackwardsPointProcess> bwdPointsInput = new Input<>("bwdPoints", "the points in the point process", Input.Validate.OPTIONAL);
    final public Input<BackwardsSchedule> bwdDisasterTimesInput = new Input<>("bwdDisasterTimes", "the times at which there was a disaster", Input.Validate.OPTIONAL);
    final public Input<BackwardsCounts> bwdDisasterCountsInput = new Input<>("bwdDisasterCounts", "the counts of removed individuals at a disaster", Input.Validate.OPTIONAL);
    final public Input<BackwardsSchedule> bwdCatastropheTimesInput = new Input<>("bwdCatastropheTimes", "the times at which there was a catastrophe", Input.Validate.OPTIONAL);

    private Tree tree;
    private BackwardsSchedule bwdCatastropheTimes;
    private BackwardsPointProcess bwdPoints;
    private BackwardsSchedule bwdDisasterTimes;
    private BackwardsCounts bwdDisasterCounts;
    private double rootLength;
    private int treeNodeCount;
    private Node[] treeNodes;
    private double[] treeNodeTimes;
    private int[] treeNodeOutdegree;
    private double[] fwdCatastropheTimes;

    private double totalTimeSpan; // the duration of time from the origin until the most recent observation.

    public TreeWithBackwardsPointProcess() {
        super();
    }

    public TreeWithBackwardsPointProcess(RealParameter rootLength,
                                         Tree tree,
                                         BackwardsPointProcess bwdPoints,
                                         BackwardsSchedule bwdDisasterTimes,
                                         BackwardsCounts bwdDisasterCounts,
                                         BackwardsSchedule bwdCatastropheTimes) {
        init(rootLength, tree, bwdPoints, bwdDisasterTimes, bwdDisasterCounts, bwdCatastropheTimes);
    }

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        treeNodes = tree.getNodesAsArray();
        treeNodeCount = tree.getNodeCount();
        treeNodeTimes = new double[treeNodeCount + 1]; // one extra time for the origin to tMCRA.
        treeNodeOutdegree = new int[treeNodeCount];
        rootLength = rootLengthInput.get().getDoubleValues()[0];
        bwdPoints = bwdPointsInput.get();
        bwdDisasterTimes = bwdDisasterTimesInput.get();
        bwdDisasterCounts = bwdDisasterCountsInput.get();
        bwdCatastropheTimes = bwdCatastropheTimesInput.get();

        calculateIntervals();
        intervalsKnown = false;
    }

    /**
     * Stolen from coalescent.TreeIntervals because it is protected there... this function mutates the times array and
     * returns nothing. Note that the times are go backwards from the present. The very last time goes all the way back
     * to the origin.
     *
     * @param rootLength
     * @param times
     * @param childCounts
     */
    protected void collectTimes(double rootLength, double[] times, int[] childCounts) {
        Node node;
        double maxHeight = 0;
        double currHeight;
        for (int i = 0; i < treeNodeCount; i++) {
            node = treeNodes[i];
            currHeight = node.getHeight();
            if (currHeight > maxHeight)
                maxHeight = currHeight;
            times[i] = currHeight;
            childCounts[i] = node.isLeaf() ? 0 : 2;
        }
        times[treeNodeCount] = maxHeight + rootLength;
    }

    /**
     * Calculate the intervals.
     */
    protected void calculateIntervals() {
        collectTimes(rootLength, treeNodeTimes, treeNodeOutdegree);
        // we need to reverse and shift the times in the tree, so it uses forward-time starting from the origin. The
        // totalTreeHeight is the length from the origin to the last tip.
        double totalTreeHeight = treeNodeTimes[treeNodeCount];
        for (int i = 0; i < treeNodeCount; i++) {
            treeNodeTimes[i] = totalTreeHeight - treeNodeTimes[i];
        }
        int[] treeJxs = new int[treeNodeCount];
        HeapSort.sort(treeNodeTimes, treeJxs);

        int[] catastropheSizes;
        if (bwdCatastropheTimes != null) {
            // we make a list of the (forward) times at which there was a catastrophe so that we can collect this
            // information out of the tree later.
            fwdCatastropheTimes = bwdCatastropheTimes.valuesInput.get().stream().mapToDouble(x -> (Double) (totalTreeHeight - x)).sorted().toArray();
            catastropheSizes = new int[fwdCatastropheTimes.length];
            Arrays.fill(catastropheSizes, 0);
            measureCatastrophes(totalTreeHeight, catastropheSizes);
        } else {
            // the catastrophe times are null which is used as a signal that there were no catastrophes.
            fwdCatastropheTimes = new double[0];
            catastropheSizes = new int[0];
        }

        // we need to handle the edge case where there are no backwards points in which case we do not expect this input
        // to be provided.
        double[] fwdPointTimes, fwdDisasterTimes;
        if (bwdPoints != null && bwdDisasterTimes != null) {
            fwdPointTimes = bwdPoints.valuesInput.get().stream().mapToDouble(x -> totalTreeHeight - x).sorted().toArray();
            fwdDisasterTimes = bwdDisasterTimes.valuesInput.get().stream().mapToDouble(x -> totalTreeHeight - x).sorted().toArray();
            double newestPointTime, newestDisasterTime,
                    oldestPointTime, oldestDisasterTime,
                    newestUnsequenced, oldestUnsequenced;
            // we have to use two streams here because the first one gets exhausted while calculating the min.
            newestPointTime = Arrays.stream(bwdPoints.getDoubleValues()).min().getAsDouble();
            newestDisasterTime = Arrays.stream(bwdDisasterTimes.getDoubleValues()).min().getAsDouble();
            newestUnsequenced = Math.min(newestPointTime, newestDisasterTime);
            oldestPointTime = Arrays.stream(bwdPoints.getDoubleValues()).max().getAsDouble();
            oldestDisasterTime = Arrays.stream(bwdDisasterTimes.getDoubleValues()).max().getAsDouble();
            oldestUnsequenced = Math.max(oldestPointTime, oldestDisasterTime);

            if (newestUnsequenced < 0 && oldestUnsequenced > totalTreeHeight) {
                setTotalTimeSpan(oldestUnsequenced - newestUnsequenced);
            } else if (newestUnsequenced >= 0 && oldestUnsequenced > totalTreeHeight) {
                setTotalTimeSpan(oldestUnsequenced);
            } else if (newestUnsequenced < 0 && oldestUnsequenced <= totalTreeHeight) {
                setTotalTimeSpan(totalTreeHeight - newestUnsequenced);
            } else { // if (newestUnsequenced >= 0 && oldestUnsequenced <= totalTreeHeight) {
                setTotalTimeSpan(totalTreeHeight);
            }
        } else if (bwdPoints != null) {
            fwdPointTimes = bwdPoints.valuesInput.get().stream().mapToDouble(x -> totalTreeHeight - x).sorted().toArray();
            fwdDisasterTimes = new double[]{};

            // we have to use two streams here because the first one gets exhausted while calculating the min.
            double newestPointTime, oldestPointTime;
            newestPointTime = Arrays.stream(bwdPoints.getDoubleValues()).min().getAsDouble();
            oldestPointTime = Arrays.stream(bwdPoints.getDoubleValues()).max().getAsDouble();

            if (newestPointTime < 0 && oldestPointTime > totalTreeHeight) {
                setTotalTimeSpan(oldestPointTime - newestPointTime);
            } else if (newestPointTime >= 0 && oldestPointTime > totalTreeHeight) {
                setTotalTimeSpan(oldestPointTime);
            } else if (newestPointTime < 0 && oldestPointTime <= totalTreeHeight) {
                setTotalTimeSpan(totalTreeHeight - newestPointTime);
            } else { // if (newestPointTime >= 0 && oldestPointTime <= totalTreeHeight) {
                setTotalTimeSpan(totalTreeHeight);
            }
        } else if (bwdDisasterTimes != null) {
            throw new RuntimeException("not implemented yet.");
        } else {
            fwdPointTimes = new double[]{};
            fwdDisasterTimes = new double[]{};
            setTotalTimeSpan(totalTreeHeight);
        }

        // The number of intervals needs to account for catastrophes where there are no sequences collected and
        // catastrophes where multiple leaves correspond to a single interval.
        intervalCount = fwdPointTimes.length +
                treeNodeCount - Arrays.stream(catastropheSizes).sum() +
                fwdCatastropheTimes.length + fwdDisasterTimes.length;
        intervals = new double[intervalCount];
        intervalTypes = new EventType[intervalCount];
        EventType birthEvent = new EventType("birth", OptionalInt.empty());
        EventType sampleEvent = new EventType("sample", OptionalInt.empty());
        EventType occurrenceEvent = new EventType("occurrence", OptionalInt.empty());

        int treeJx = 0;
        int pointJx = 0;
        int catastJx = 0;
        int disastJx = 0;
        // while looping it is helpful to have a reference to the tree event time, point event time, etc.
        double treeET, pointET, catastET, disastET;
        double currTime = 0.0; // TODO we should not really be hardcoding that the origin is zero....
        int intIx = 0;
        while (intIx < intervalCount) {
            // if the tree events have been exhausted set the next event time to positive
            // infinity so that it cannot appear as the next observation. Do the same the point process and the
            // disasters.
            treeET = (treeJx < treeJxs.length) ? treeNodeTimes[treeJxs[treeJx]] : Double.POSITIVE_INFINITY;
            pointET = (pointJx < fwdPointTimes.length) ? fwdPointTimes[pointJx] : Double.POSITIVE_INFINITY;
            // we keep track of the time of the next catastrophe so that it is easier to check how many sequences where
            // sampled and when it occurred.
            catastET = (catastJx < fwdCatastropheTimes.length) ? fwdCatastropheTimes[catastJx] : Double.POSITIVE_INFINITY;
            disastET = (disastJx < fwdDisasterTimes.length) ? fwdDisasterTimes[disastJx] :  Double.POSITIVE_INFINITY;

            if ((treeET == pointET && Double.isFinite(treeET))
                || (treeET == disastET && Double.isFinite(treeET))
                || (pointET == disastET && Double.isFinite(pointET))) {
                throw new IllegalArgumentException(
                        "It appears that a sequenced and unsequenced event occurred at the same time.");
            }

            if (treeET < pointET && treeET < disastET) {
                intervals[intIx] = treeET - currTime;
                currTime = treeET;
                if (Math.abs(treeET - catastET) > 1e-9) {
                    if (treeNodeOutdegree[treeJxs[treeJx]] == 2) {
                        intervalTypes[intIx] = birthEvent;
                    } else if (treeNodeOutdegree[treeJxs[treeJx]] == 0) {
                        intervalTypes[intIx] = sampleEvent;
                    } else {
                        throw new IllegalArgumentException("Non-binary tree provided to TreeWithPointProcess");
                    }
                    treeJx++;
                } else {
                    intervalTypes[intIx] = new EventType("catastrophe", OptionalInt.of(catastropheSizes[catastJx]));
                    catastJx++;
                    while (Math.abs(treeET - catastET) < 1e-9) {
                        treeJx++;
                        treeET = (treeJx < treeJxs.length) ? treeNodeTimes[treeJxs[treeJx]] : Double.POSITIVE_INFINITY;
                    }
                }
            } else if (pointET < disastET) {
                // the next event must have been an occurrence.
                intervals[intIx] = pointET - currTime;
                currTime = pointET;
                intervalTypes[intIx] = occurrenceEvent;
                pointJx++;
            } else {
                // the next event must have been a disaster.
                intervals[intIx] = disastET - currTime;
                currTime = disastET;
                intervalTypes[intIx] = new EventType("disaster", OptionalInt.of(bwdDisasterCounts.getNativeValue(fwdDisasterTimes.length - disastJx - 1)));
                disastJx++;
            }
            intIx++;
        }

        intervalsKnown = true;
    }

    /**
     * Look at the tree and detect catastrophes and mutate the array of sizes to contain these counts.
     *
     * @param maxTime
     * @param catastropheSizes
     * @see TreeWithBackwardsPointProcess#collectTimes
     */
    private void measureCatastrophes(double maxTime, int[] catastropheSizes) {
        double nodeTime;

        for (Node node : treeNodes) {
            if (node.isLeaf()) {
                nodeTime = maxTime - node.getHeight();
                for (int j = 0; j < this.fwdCatastropheTimes.length; j++) {
                    if (Math.abs(this.fwdCatastropheTimes[j] - nodeTime) < 1e-8) {
                        catastropheSizes[j] += 1;
                    }
                }
            }
        }

        if (Arrays.stream(catastropheSizes).sum() == 0) {
            StringBuilder nodeTimesSB = new StringBuilder();
            for (Node node : treeNodes) {
                if (node.isLeaf()) {
                    nodeTimesSB.append("\n\t\t").append(maxTime - node.getHeight());
                }
            }
            throw new RuntimeException(
                    "\nIt appears that no samples were obtained in any of the catastrophies, please " +
                            "\ndouble check that the timing of these has been specified correctly in the XML. " +
                            "\nFor debugging purposes here are some variables:" +
                            "\n\tThe current catastrophe times are " + this.fwdCatastropheTimes +
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
        if (i < 0 || i >= intervalCount)
            throw new IllegalArgumentException("intervalCount is " + intervalCount + " but interval " + i + " requested.");
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

    /**
     * This is the total amount of time from the origin of the process until the time of the final observation.
     *
     * @return total duration of time from origin to last observation
     */
    public double getTotalTimeSpan() {
        return totalTimeSpan;
    }

    public void setTotalTimeSpan(double totalTimeSpan) {
        this.totalTimeSpan = totalTimeSpan;
    }

    /**
     * Predicate to check that the values are non-decreasing.
     */
    static boolean isMonotonicIncreasing(double[] xs) {
        if (xs.length > 1) {
            for (int i = 1; i < xs.length; i++) {
                if (xs[i-1] > xs[i]) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * Predicate to check that the values are non-increasing.
     */
    static boolean isMonotonicDecreasing(double[] xs) {
        if (xs.length > 1) {
            for (int i = 1; i < xs.length; i++) {
                if (xs[i-1] < xs[i]) {
                    return false;
                }
            }
        }
        return true;
    }
}