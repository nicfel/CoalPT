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
        " the framework of Mueller (2022).")
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
//    	System.out.println("====");
//    	System.out.println(intervals.networkInput.get());
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
					logP += plasmidTransfer(event);
					break;
			}

       		if (logP==Double.NEGATIVE_INFINITY)
       			break;

        	prevEvent = event;
        }        
//    	System.out.println(logP + " " + storedLogP);
		return logP;
    }
    
	private double plasmidTransfer(NetworkEvent event) {
		if (intervals.hasMultipleRates) {
			if (event.segsLeft.cardinality()==0 || event.segsRight.cardinality()==0) {
				return Double.NEGATIVE_INFINITY;
			}
			
			if (event.segsLeft.cardinality()>1) {
				return Math.log(intervals.plasmidTransferRate.getArrayValue(event.segsRight.nextSetBit(0)-1));
			}else if (event.segsRight.cardinality()>1){
				return Math.log(intervals.plasmidTransferRate.getArrayValue(event.segsLeft.nextSetBit(0)-1));
				
			}else if (event.segsLeft.get(0)){
				return Math.log(intervals.plasmidTransferRate.getArrayValue(event.segsRight.nextSetBit(0)-1));
			}else if (event.segsRight.get(0)){
				return Math.log(intervals.plasmidTransferRate.getArrayValue(event.segsLeft.nextSetBit(0)-1));
			}else {
				double rate = intervals.plasmidTransferRate.getArrayValue(event.segsLeft.nextSetBit(0)-1) +
						intervals.plasmidTransferRate.getArrayValue(event.segsRight.nextSetBit(0)-1);
				return Math.log(rate/2);
			}

		}else {
			return Math.log(intervals.plasmidTransferRate.getArrayValue());
		}
	}

	private double coalesce(NetworkEvent event) {

		return Math.log(1.0/populationFunction.getPopSize(event.time));
	}

	private double intervalContribution(NetworkEvent prevEvent, NetworkEvent nextEvent) {

        double result = 0.0;

        result += -prevEvent.totalReassortmentObsProb
                * (nextEvent.time - prevEvent.time);

		result += -0.5*prevEvent.lineages*(prevEvent.lineages-1)
                * populationFunction.getIntegral(prevEvent.time, nextEvent.time);
		
		return result;
	}
	
//    @Override
//    protected boolean requiresRecalculation() {    	   	
//    	if (((CalculationNode) populationFunction).isDirtyCalculation())
//    		return true;
//    	
//        return super.requiresRecalculation();
//    }
    

}
