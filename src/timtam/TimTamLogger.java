package timtam;


import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;

import java.io.PrintStream;

public class TimTamLogger extends BEASTObject implements Loggable {

    // TODO This should really include an additional input to specify the parameterisation of the prevelence estimate to return: r and p or the mean and variance.
    final public Input<TimTam> timTamInput = new Input<>("timtam", "the prevalence estimate generated while calculating the likelihood using timtam.", Input.Validate.REQUIRED);

    @Override
    public void init(PrintStream out) {
        final TimTam tt = timTamInput.get();
        out.print(tt.getID() + ".prevalence.lnR\t");
        out.print(tt.getID() + ".prevalence.lnP\t");
        out.print(tt.getID() + ".prevalence.lnMean\t");
        out.print(tt.getID() + ".prevalence.lnVariance\t");
    }

    @Override
    public void log(long sample, PrintStream out) {
        final TimTam tt = timTamInput.get();
        out.print(tt.getTimTamNegBinom().getLnR() + "\t");
        out.print(tt.getTimTamNegBinom().getLnP() + "\t");
        out.print(tt.getTimTamNegBinom().getLnMean() + "\t");
        out.print(tt.getTimTamNegBinom().getLnVariance() + "\t");
    }

    @Override
    public void close(PrintStream out) {

    }

    @Override
    public void initAndValidate() {

    }
}
