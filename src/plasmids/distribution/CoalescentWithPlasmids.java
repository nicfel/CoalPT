package plasmids.distribution;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.coalescent.PopulationFunction;
import coalre.distribution.NetworkEvent;

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
	
	public Input<Boolean> conditionOnCoalescentEventsInput = new Input<>(
	        "conditionOnCoalescentEvents",
            "if true, only coalescent events are allowed after the .",
            true);	
	
	public Input<Double> maxHeightRatioInput = new Input<>(
	        "maxHeightRatio",
            "if specified, above the ratio, only coalescent events are allowed.", 1.25);	

	public Input<Double> redFactorInput = new Input<>(
	        "redFactor",
            "by how much the recombination rate should be reduced after reaching the maxHeightRatio.", 0.0);	

	

    public PopulationFunction populationFunction;

    public PlasmidNetworkIntervals intervals;

    @Override
    public void initAndValidate(){
        populationFunction = populationFunctionInput.get();
        intervals = networkIntervalsInput.get();
    }

//    int iii=0;
    public double calculateLogP() {
    	logP = 0;
    	
    	// Calculate tree intervals
    	List<NetworkEvent> networkEventList = intervals.getNetworkEventList();

    	NetworkEvent prevEvent = null;
    	
    	// get the mrca of all loci trees
    	double lociMRCA = conditionOnCoalescentEventsInput.get() ? intervals.getMaxSegmentTreeHeight()*maxHeightRatioInput.get() : Double.POSITIVE_INFINITY;
    	

    	for (NetworkEvent event : networkEventList) {
        	if (prevEvent != null)
        		logP += intervalContribution(prevEvent, event, lociMRCA);

        	switch (event.type) {
				case COALESCENCE:
					logP += coalesce(event);
					break;

				case SAMPLE:
					break;

				case REASSORTMENT:
					logP += plasmidTransfer(event, lociMRCA);
					break;
			}

       		if (logP==Double.NEGATIVE_INFINITY)
       			break;

        	prevEvent = event;
        }        
		return logP;
    }
    
	private double plasmidTransfer(NetworkEvent event, double lociMRCA) {
		if (intervals.hasMultipleRates) {
			if (event.segsLeft.cardinality()==0 || event.segsRight.cardinality()==0) {
				return Double.NEGATIVE_INFINITY;
			}
			
			if (event.time>lociMRCA) {
			
				if (event.segsLeft.cardinality()>1) {
					return Math.log(intervals.plasmidTransferRate.getArrayValue(event.segsRight.nextSetBit(0)-1) * redFactorInput.get() );
				}else if (event.segsRight.cardinality()>1){
					return Math.log(intervals.plasmidTransferRate.getArrayValue(event.segsLeft.nextSetBit(0)-1) * redFactorInput.get());
					
				}else if (event.segsLeft.get(0)){
					return Math.log(intervals.plasmidTransferRate.getArrayValue(event.segsRight.nextSetBit(0)-1) * redFactorInput.get());
				}else if (event.segsRight.get(0)){
					return Math.log(intervals.plasmidTransferRate.getArrayValue(event.segsLeft.nextSetBit(0)-1) * redFactorInput.get());
				}else {
					double rate = intervals.plasmidTransferRate.getArrayValue(event.segsLeft.nextSetBit(0)-1) +
							intervals.plasmidTransferRate.getArrayValue(event.segsRight.nextSetBit(0)-1);
					
					return Math.log(rate/2 * redFactorInput.get());
				}
			}else {
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

			}

		}else {
			if (event.time>lociMRCA) {
				return Math.log(intervals.plasmidTransferRate.getArrayValue()*redFactorInput.get());
			}else {
				return Math.log(intervals.plasmidTransferRate.getArrayValue());
			}
		}
	}

	private double coalesce(NetworkEvent event) {

		return Math.log(1.0/populationFunction.getPopSize(event.time));
	}

	private double intervalContribution(NetworkEvent prevEvent, NetworkEvent nextEvent, double lociMRCA) {

        double result = 0.0;
        
        double scale_factor = 1.0;
        
        if (nextEvent.time>lociMRCA) {   
        	if (prevEvent.time<lociMRCA) 
            	scale_factor = ((lociMRCA-prevEvent.time) +  redFactorInput.get() * (nextEvent.time-lociMRCA))/(nextEvent.time-prevEvent.time);
        	else
        		scale_factor = redFactorInput.get();
        }	        

        

        result += -prevEvent.totalReassortmentObsProb
                * (nextEvent.time - prevEvent.time) * scale_factor;

		result += -0.5*prevEvent.lineages*(prevEvent.lineages-1)
                * populationFunction.getIntegral(prevEvent.time, nextEvent.time);
		
		return result;
	}
	
    @Override
    protected boolean requiresRecalculation() {    	   	
    	if (((CalculationNode) populationFunction).isDirtyCalculation())
    		return true;
    	
        return super.requiresRecalculation();
    }	
}
