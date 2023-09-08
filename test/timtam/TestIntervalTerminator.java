package timtam;

import org.junit.Test;

import java.util.Arrays;
import java.util.Objects;
import java.util.OptionalInt;

import static org.junit.Assert.assertTrue;

public class TestIntervalTerminator {

    @Test
    public void testOrdering() {

        IntervalTerminator it1 =
                new IntervalTerminator("birth", 1.0);
        IntervalTerminator it2 =
                new IntervalTerminator("birth", 0.5);
        IntervalTerminator it3 =
                new IntervalTerminator("birth", 0.25);

        IntervalTerminator[] intTerms = {it3, it1, it2};
        Arrays.sort(intTerms);

        assertTrue(Objects.equals(intTerms[0].getBwdTime(), it1.getBwdTime()));
        assertTrue(Objects.equals(intTerms[1].getBwdTime(), it2.getBwdTime()));
        assertTrue(Objects.equals(intTerms[2].getBwdTime(), it3.getBwdTime()));
    }

    @Test(expected = RuntimeException.class)
    public void testGuardAgainstNonsense() {
        IntervalTerminator it =
                new IntervalTerminator("occurrence", 1.0, 2);
    }
}
