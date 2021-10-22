package beast.evolution.tree.birthdeath;

/**
 * Specifies that ways that an interval in a reconstructed tree can end.
 */
public enum EventType {

    /**
     * Denotes a bifuration of the in the reconstructed tree.
     */
    BIRTH("birth"),

    /**
     * Denotes a leaf of the reconstructed tree.
     */
    SAMPLE("sample"),

    /**
     * Denotes an unsequenced observation.
     */
    OCCURRENCE("occurrence");

    EventType(String name) { this.name = name; }

    @Override
    public String toString() { return name; }

    private final String name;
}
