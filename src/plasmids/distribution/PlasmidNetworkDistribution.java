package plasmids.distribution;

import beast.core.CalculationNode;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;

import java.util.List;
import java.util.Random;

public class PlasmidNetworkDistribution extends Distribution {
    public Input<PlasmidNetworkIntervals> networkIntervalsInput = new Input<>("networkIntervals",
            "Structured Intervals for a phylogenetic beast tree", Validate.REQUIRED);

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
    }
    

}
