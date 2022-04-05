package timtam;

import java.util.Objects;
import java.util.OptionalInt;
import java.util.Set;

public class TimTamIntervalTerminator implements Comparable<TimTamIntervalTerminator> {

    static final Set<String> validTerminatorNames =
            Set.of("birth", "sample", "occurrence", "catastrophe", "disaster", "rateChange");

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

    TimTamIntervalTerminator(String type, Double bwdTime, OptionalInt count) {
        if (validTerminatorNames.contains(type)) {
            this.type = type;
            this.bwdTime = bwdTime;
            if (count.isPresent() & (Objects.equals(type, "catastrophe") | Objects.equals(type, "disaster"))) {
                this.count = count.getAsInt();
            } else if (count.isEmpty() & (!Objects.equals(type, "catastrophe")) & (!Objects.equals(type, "disaster"))) {
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
