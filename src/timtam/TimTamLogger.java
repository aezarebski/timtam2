package timtam;


import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;

import java.io.PrintStream;

public class TimTamLogger extends BEASTObject implements Loggable {

    final public Input<TimTam> timTamInput = new Input<TimTam>("timtam", "the prevalence distribution generated while calculating the likelihood using timtam.", Input.Validate.REQUIRED);

    @Override
    public void init(PrintStream out) {
        final TimTam tt = timTamInput.get();
        out.print(tt.getID() + ".prevalence.mean\t");
        out.print(tt.getID() + ".prevalence.variance\t");
    }

    @Override
    public void log(long sample, PrintStream out) {
        final TimTam tt = timTamInput.get();
        out.print(Math.exp(tt.getTimTamNegBinom().getLnMean()) + "\t");
        out.print(Math.exp(tt.getTimTamNegBinom().getLnVariance()) + "\t");
    }

    @Override
    public void close(PrintStream out) {

    }

    @Override
    public void initAndValidate() {

    }
}
