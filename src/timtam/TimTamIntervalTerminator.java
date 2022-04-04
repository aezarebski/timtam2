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
            return count.getAsInt();
        } else {
            throw new RuntimeException("TimTamIntervalTerminator, " + this.type + " does not have a count.");
        }
    }

    public Double getBwdTime() {
        return bwdTime;
    }

    private final String type;
    private final OptionalInt count;
    private final Double bwdTime;

    TimTamIntervalTerminator(String type, Double bwdTime, OptionalInt count) {
        if (validTerminatorNames.contains(type)) {
            this.type = type;
            this.bwdTime = bwdTime;
            this.count = count;
        } else {
            throw new IllegalArgumentException("Unexpected interval terminator type: " + type);
        }
    }

    @Override
    public int compareTo(TimTamIntervalTerminator o) {
        return o.bwdTime.compareTo(this.bwdTime);
    }
}
