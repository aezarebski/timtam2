package timtam;

import java.util.*;

public class IntervalTerminator implements Comparable<IntervalTerminator> {

    private String type;
    private int count;
    private boolean shouldHaveCount;
    private Double bwdTime;
    private final Set<String> validTerminatorNames =
            new HashSet<>(Arrays.asList(
                    "birth",
                    "sample",
                    "occurrence",
                    "catastrophe",
                    "disaster",
                    "paramValueChange",
                    "historyEstimate"
            ));

    public IntervalTerminator() {
    }

    /**
     * Construct a description of why an interval terminated. This is used by {@link TimTam}.
     *
     * @param type a string indicating the reason the interval ended. This must be one of the values in {@link IntervalTerminator#validTerminatorNames}.
     * @param bwdTime the time the interval ended measured backwards from the time of the last sequenced observation.
     * @param count the size of the event (typically the number of lineages involved in the event).
     *
     * <b>NOTE:</b> If there is no count given to this constructor will set it to a dummy value of -1.
     *
     */
    public IntervalTerminator(String type, Double bwdTime, OptionalInt count) {
        if (validTerminatorNames.contains(type)) {
            this.type = type;
            this.bwdTime = bwdTime;
            this.shouldHaveCount =
                    Objects.equals(type, "catastrophe") |
                    Objects.equals(type, "disaster") |
                    Objects.equals(type, "historyEstimate");
            if (count.isPresent() & this.shouldHaveCount) {
                this.count = count.getAsInt();
            } else if ((!count.isPresent()) & (!this.shouldHaveCount)) {
                this.count = -1;
            } else {
                throw new RuntimeException(
                        "count attribute cannot be set to "+ count + " in TimTamIntervalTerminator if type is " + type
                );
            }
        } else {
            throw new IllegalArgumentException("Unexpected interval terminator type: " + type);
        }
    }

    public String getType() {
        return type;
    }

    public int getCount() {
        if (this.shouldHaveCount) {
            return this.count;
        } else {
            throw new RuntimeException("TimTamIntervalTerminator, " + this.type + " does not have a count.");
        }
    }

    public Double getBwdTime() {
        return bwdTime;
    }

    /**
     * Update the data associated with this object.
     *
     * @param type a string indicating the reason the interval ended. This must be one of the values in {@link IntervalTerminator#validTerminatorNames}.
     * @param bwdTime the time the interval ended measured backwards from the time of the last sequenced observation.
     * @param count the size of the event (typically the number of lineages involved in the event).
     *
     * <b>NOTE:</b> If there is no count given to this constructor will set it to a dummy value of -1.
     */
    public void setTypeTimeAndCount(String type, Double bwdTime, OptionalInt count) {
        if (validTerminatorNames.contains(type)) {
            this.type = type;
            this.bwdTime = bwdTime;
            this.shouldHaveCount =
                    Objects.equals(type, "catastrophe") |
                    Objects.equals(type, "disaster") |
                    Objects.equals(type, "historyEstimate");
            if (count.isPresent() & this.shouldHaveCount) {
                this.count = count.getAsInt();
                if (this.count < 0) {
                    throw new RuntimeException("count cannot be "+ count + " in TimTamIntervalTerminator if type is " + type);
                }
            } else if ((!count.isPresent()) & (!this.shouldHaveCount)) {
                this.count = -1;
            } else {
                throw new RuntimeException("count cannot be "+ count + " in TimTamIntervalTerminator if type is " + type);
            }
        } else {
            throw new IllegalArgumentException("Unexpected interval terminator type: " + type);
        }
    }

    @Override
    public int compareTo(IntervalTerminator o) {
        return o.bwdTime.compareTo(this.bwdTime);
    }
}
