package timtam;


import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;

import java.io.PrintStream;

/**
 * <p>Class to log the distribution of the prevalence as computed by {@link
 * TimTam}.</p>
 *
 * @author Alexander E. Zarebski
 */
public class TimTamLogger extends BEASTObject implements Loggable {

    final public Input<TimTam> timTamInput = new Input<>("timtam", "the prevalence distribution generated while calculating the likelihood using timtam.", Input.Validate.REQUIRED);

    final public Input<Boolean> reportFirstIntervalInput = new Input<>("reportFirstInterval", "if is true log the length of the first interval. The default value is false.", false);

    private boolean reportFirstInterval;
    private TimTam timTam;

    @Override
    public void init(PrintStream out) {
        out.print(this.timTam.getID() + ".prevalence.mean\t");
        out.print(this.timTam.getID() + ".prevalence.variance\t");
        if (reportFirstInterval) {
            out.println(this.timTam.getID() + "first.interval\t");
        }
    }

    @Override
    public void log(long sample, PrintStream out) {
        out.print(Math.exp(this.timTam.getTimTamNegBinom().getLnMean()) + "\t");
        out.print(Math.exp(this.timTam.getTimTamNegBinom().getLnVariance()) + "\t");
        if (reportFirstInterval) {
            out.print(this.timTam.getFirstIntervalDuration() + "\t");
        }
    }

    @Override
    public void close(PrintStream out) {

    }

    @Override
    public void initAndValidate() {
        this.reportFirstInterval = reportFirstIntervalInput.get();
        this.timTam = timTamInput.get();
    }
}
