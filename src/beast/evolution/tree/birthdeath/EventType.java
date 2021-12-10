package beast.evolution.tree.birthdeath;

import java.util.Arrays;
import java.util.List;
import java.util.OptionalInt;

/**
 * Specifies that ways that an interval of time without an observation can end.
 *
 * <p>This uses a class rather than an enumeration because we do not know ahead
 * of time how many samples may have been generated in a catastrophe (which is a
 * scheduled sequenced sample) or a dissaster (which is a scheduled unsequenced
 * sample)</p>
 */
public class EventType {

    static final List<String> eventTypes = Arrays.asList("birth", "sample", "occurrence", "catastrophe", "disaster");

    private final String type;

    private final OptionalInt count;

    EventType(String type, OptionalInt count) {
        if (eventTypes.contains(type)) {
            this.type = type;
            this.count = count;
        } else {
            throw new IllegalArgumentException("Unexpected event type: " + type);
        }
    }

    @Override
    public String toString() { return type; }

    public OptionalInt getCount() {
        return count;
    }
}
