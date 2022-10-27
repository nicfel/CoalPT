package plasmids.ratemodel;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.coalescent.PopulationFunction;

public class EarlyLateRate extends BranchRateModel.Base {

	final public Input<RealParameter> lateRateInput = new Input<>("lateRate",
			"the Clock rates through time.", Input.Validate.REQUIRED);
	
	final public Input<Double> finalTimeInput = new Input<>("finalTime",
			"the time when to switch to the late clock rate.", Input.Validate.REQUIRED);


	@Override
	public void initAndValidate() {
	}
	

	@Override
	public double getRateForBranch(Node node) {			
		if (node.isRoot()) {
			return 0.0;		
		}	
		if (node.getParent().getHeight()<finalTimeInput.get()) {
			return meanRateInput.get().getArrayValue();
		}else if (node.getHeight()>finalTimeInput.get()) {
			return lateRateInput.get().getArrayValue();
		}else {
			double part1 = finalTimeInput.get()-node.getHeight();
			double part2 = node.getParent().getHeight()-finalTimeInput.get();
			return (part1*meanRateInput.get().getArrayValue() +part2*lateRateInput.get().getArrayValue())/(part1+part2);
		}

	}


	@Override
	protected boolean requiresRecalculation() {
		if (((CalculationNode) lateRateInput.get()).isDirtyCalculation()) {
			// this is only called if any of its inputs is dirty, hence we need to recompute
			return true;
		}
		return super.requiresRecalculation();
	}


}
