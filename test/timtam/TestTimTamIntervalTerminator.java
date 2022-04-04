package timtam;

import org.junit.Test;

import java.util.Arrays;
import java.util.Objects;
import java.util.OptionalInt;

import static org.junit.Assert.assertTrue;

public class TestTimTamIntervalTerminator {

    @Test
    public void testOrdering() {

        TimTamIntervalTerminator it1 =
                new TimTamIntervalTerminator("birth", 1.0, OptionalInt.empty());
        TimTamIntervalTerminator it2 =
                new TimTamIntervalTerminator("birth", 0.5, OptionalInt.empty());
        TimTamIntervalTerminator it3 =
                new TimTamIntervalTerminator("birth", 0.25, OptionalInt.empty());

        TimTamIntervalTerminator[] intTerms = {it3, it1, it2};
        Arrays.sort(intTerms);

        assertTrue(Objects.equals(intTerms[0].getBwdTime(), it1.getBwdTime()));
        assertTrue(Objects.equals(intTerms[1].getBwdTime(), it2.getBwdTime()));
        assertTrue(Objects.equals(intTerms[2].getBwdTime(), it3.getBwdTime()));
    }
}
