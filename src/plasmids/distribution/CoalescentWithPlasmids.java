package plasmids.distribution;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.evolution.tree.coalescent.PopulationFunction;
import coalre.distribution.NetworkDistribution;
import coalre.distribution.NetworkEvent;
import coalre.distribution.NetworkIntervals;

import java.util.List;


/**
 * @author Nicola Felix Mueller
 */

@Description("Calculates the probability of a core genomes + plasmids network under" +
        " the framework of Mueller (2021).")
public class CoalescentWithPlasmids extends PlasmidNetworkDistribution {

	public Input<PopulationFunction> populationFunctionInput = new Input<>(
	        "populationModel",
            "Population model.",
            Input.Validate.REQUIRED);
	

    public PopulationFunction populationFunction;

    public PlasmidNetworkIntervals intervals;

    @Override
    public void initAndValidate(){
        populationFunction = populationFunctionInput.get();
        intervals = networkIntervalsInput.get();
    }

    public double calculateLogP() {
    	logP = 0;
    	// Calculate tree intervals
    	List<NetworkEvent> networkEventList = intervals.getNetworkEventList();

    	NetworkEvent prevEvent = null;

    	for (NetworkEvent event : networkEventList) {
        	if (prevEvent != null)
        		logP += intervalContribution(prevEvent, event);

        	switch (event.type) {
				case COALESCENCE:
					logP += coalesce(event);
					break;

				case SAMPLE:
					break;

				case REASSORTMENT:
					logP += reassortment(event);
					break;
			}

       		if (logP==Double.NEGATIVE_INFINITY)
       			break;

        	prevEvent = event;
        }        
		return logP;
    }
    
	private double reassortment(NetworkEvent event) {
		
		double binomval = Math.pow(intervals.getBinomialProb(), event.segsSortedLeft)
				* Math.pow(1-intervals.getBinomialProb(), event.segsToSort-event.segsSortedLeft) 
				+ Math.pow(intervals.getBinomialProb(), event.segsToSort-event.segsSortedLeft)
				* Math.pow(1-intervals.getBinomialProb(), event.segsSortedLeft); 
				        
		return Math.log(intervals.plasmidTransferRate.getArrayValue())
				+ Math.log(binomval);
	}

	private double coalesce(NetworkEvent event) {

		return Math.log(1.0/populationFunction.getPopSize(event.time));
	}

	private double intervalContribution(NetworkEvent prevEvent, NetworkEvent nextEvent) {

        double result = 0.0;

        result += -intervals.plasmidTransferRate.getArrayValue() * prevEvent.totalReassortmentObsProb
                * (nextEvent.time - prevEvent.time);

		result += -0.5*prevEvent.lineages*(prevEvent.lineages-1)
                * populationFunction.getIntegral(prevEvent.time, nextEvent.time);
		
		return result;
	}
	
    @Override
    protected boolean requiresRecalculation() {    	
    	if (((CalculationNode) plasmidTransferRate).isDirtyCalculation())
    		return true;
    	
    	if (((CalculationNode) populationFunction).isDirtyCalculation())
    		return true;
    	
        return super.requiresRecalculation();
    }
    

}
