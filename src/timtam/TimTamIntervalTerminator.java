package timtam;

import java.util.*;

public class TimTamIntervalTerminator implements Comparable<TimTamIntervalTerminator> {

    public String getType() {
        return type;
    }

    public int getCount() {
        if (Objects.equals(this.type, "catastrophe") | Objects.equals(this.type, "disaster")) {
            return count;
        } else {
            throw new RuntimeException("TimTamIntervalTerminator, " + this.type + " does not have a count.");
        }
    }

    public Double getBwdTime() {
        return bwdTime;
    }

    private final String type;
    private final int count;
    private final Double bwdTime;

    /**
     * Construct a description of why an interval terminated.
     *
     * @param type a string indicating the reason the interval ended.
     * @param bwdTime the time that the interval ended
     * @param count the size of the event
     */
    TimTamIntervalTerminator(String type, Double bwdTime, OptionalInt count) {
        final Set<String> validTerminatorNames =
                new HashSet<>(Arrays.asList("birth", "sample", "occurrence", "catastrophe", "disaster", "rateChange"));

        if (validTerminatorNames.contains(type)) {
            this.type = type;
            this.bwdTime = bwdTime;
            if (count.isPresent() & (Objects.equals(type, "catastrophe") | Objects.equals(type, "disaster"))) {
                this.count = count.getAsInt();
            } else if ((!count.isPresent()) & (!Objects.equals(type, "catastrophe")) & (!Objects.equals(type, "disaster"))) {
                this.count = -1;
            } else {
                throw new RuntimeException("count cannot be "+ count + " in TimTamIntervalTerminator if type is " + type);
            }
        } else {
            throw new IllegalArgumentException("Unexpected interval terminator type: " + type);
        }
    }

    @Override
    public int compareTo(TimTamIntervalTerminator o) {
        return o.bwdTime.compareTo(this.bwdTime);
    }
}
