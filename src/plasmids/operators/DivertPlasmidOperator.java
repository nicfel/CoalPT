package plasmids.operators;

import beast.util.Randomizer;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;
import coalre.operators.EmptyEdgesNetworkOperator;

import java.util.BitSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public class DivertPlasmidOperator extends EmptyEdgesNetworkOperator {

	// probability of a plasmid going left or right
	double prob;
	
    @Override
    public void initAndValidate() {
        super.initAndValidate();
    	prob = 1.0/(network.getSegmentCount()-1);
   	}

	
    @Override
    public double networkProposal() {
		throw new IllegalArgumentException("not an operator");
    }


    /**
     * Remove segments from this edge and ancestors.
     *
     * @param edge edge at which to start removal
     * @param segsToRemove segments to remove from edge and ancestors
     * @return log probability of reverse operation
     */
    double removeSegmentsFromAncestors(NetworkEdge edge, BitSet segsToRemove) {
        double logP = 0.0;

        segsToRemove = (BitSet)segsToRemove.clone();
        segsToRemove.and(edge.hasSegments);

        if (segsToRemove.isEmpty())
            return logP;

        edge.hasSegments.andNot(segsToRemove);

        if (edge.isRootEdge())
            return logP;

        if (edge.parentNode.isReassortment()) {
        	
            BitSet segsToAddLeft = (BitSet) edge.parentNode.getParentEdges().get(0).hasSegments.clone();
            BitSet segsToAddRight = (BitSet) edge.parentNode.getParentEdges().get(1).hasSegments.clone();
            
            segsToAddLeft.andNot(segsToRemove);
            segsToAddRight.andNot(segsToRemove);
            
            for (int segIdx=segsToRemove.nextSetBit(0); segIdx != -1;
                    segIdx=segsToRemove.nextSetBit(segIdx+1)) {
            	
            	if (segsToAddLeft.cardinality()==0 && segsToAddRight.cardinality()==0) { // randomly pick a path when both are empty
                    if (edge.parentNode.getParentEdges().get(0).hasSegments.get(segIdx))
                        segsToAddLeft.set(segIdx);
                    else
                        segsToAddRight.set(segIdx);

                    logP += Math.log(0.5);
            	}else if (segsToAddLeft.get(0) && segsToAddRight.cardinality()==1) { // has to follow the core genome when there is already a plasmid transferred
            		segsToAddLeft.set(segIdx);
            	}else if (segsToAddRight.get(0) && segsToAddLeft.cardinality()==1) { // has to follow the core genome when there is already a plasmid transferred
                    segsToAddRight.set(segIdx);
            	}else if (segsToAddLeft.cardinality()>1 && segsToAddRight.cardinality()==1) { // has to follow the path with more plasmids already there as only one plasmid is transferred at a time
            		segsToAddLeft.set(segIdx);
            	}else if (segsToAddRight.cardinality()>1 && segsToAddLeft.cardinality()==1) { // has to follow the path with more plasmids already there as only one plasmid is transferred at a time
            		segsToAddRight.set(segIdx);

            	}else if (segsToAddLeft.cardinality()==1 && segsToAddRight.cardinality()==1) { // as neither path has the core genome, just randomly pick a path
                    if (edge.parentNode.getParentEdges().get(0).hasSegments.get(segIdx))
                        segsToAddLeft.set(segIdx);
                    else
                        segsToAddRight.set(segIdx);
                    logP += Math.log(0.5);
            	}else if (segsToAddLeft.cardinality()==1 && segsToAddRight.cardinality()==0) { // left path has one, right has 0, doesn't matter where core goes
                    if (edge.parentNode.getParentEdges().get(1).hasSegments.get(segIdx)) {
                    	segsToAddRight.set(segIdx);
                        logP += Math.log(prob);
                    }else {
                    	segsToAddLeft.set(segIdx);
                        logP += Math.log(1-prob);
                    }
            	}else if (segsToAddRight.cardinality()==1 && segsToAddLeft.cardinality()==0) { // right path has one, left has 0, doesn't matter where core goes
                    if (edge.parentNode.getParentEdges().get(0).hasSegments.get(segIdx)) {
                    	segsToAddLeft.set(segIdx);
                        logP += Math.log(prob);
                    }else {
                    	segsToAddRight.set(segIdx);
                        logP += Math.log(1-prob);
                    }
            	}else if (segsToAddLeft.cardinality()>1 && segsToAddRight.cardinality()==0) { // left path has 2+, right has 0, core goes left, otherwise random
            		if (segIdx==0) {
                    	segsToAddLeft.set(segIdx);
            		}else if (edge.parentNode.getParentEdges().get(1).hasSegments.get(segIdx)) {
            			segsToAddRight.set(segIdx);
                        logP += Math.log(1-prob);
                    }else {
                    	segsToAddLeft.set(segIdx);
                        logP += Math.log(1-prob);
                    }
            	} else if (segsToAddRight.cardinality()>1 && segsToAddLeft.cardinality()==0) { // right path has 2+, left has 0, core goes right, otherwise random
            		if (segIdx==0) {
            			segsToAddRight.set(segIdx);
            		}else if (edge.parentNode.getParentEdges().get(0).hasSegments.get(segIdx)) {
            			segsToAddLeft.set(segIdx);
                        logP += Math.log(prob);
                    }else {
                    	segsToAddRight.set(segIdx);
                        logP += Math.log(1-prob);
                    }

            	} else {
            		throw new IllegalArgumentException("scenario unknown, should not happen");
            	}
            	
            }

            // sanity check 
            

            logP += removeSegmentsFromAncestors(edge.parentNode.getParentEdges().get(0), segsToRemove);
            logP += removeSegmentsFromAncestors(edge.parentNode.getParentEdges().get(1), segsToRemove);

        } else {

            segsToRemove.andNot(getSisterEdge(edge).hasSegments);
            logP += removeSegmentsFromAncestors(edge.parentNode.getParentEdges().get(0), segsToRemove);

        }
        
        
        return logP;
    }

    /**
     * Add segments to this edge and ancestors, while only allowing for one plasmid to be transferred at a time.
     *
     * @param edge edge at which to start addition
     * @param segsToAdd segments to add to the edge and ancestors
     * @return log probability of operation
     */
    double addSegmentsToAncestors(NetworkEdge edge, BitSet segsToAdd) {
        double logP = 0.0;

        segsToAdd = (BitSet) segsToAdd.clone();
        segsToAdd.andNot(edge.hasSegments);

        if (segsToAdd.isEmpty())
            return logP;

        edge.hasSegments.or(segsToAdd);

        if (edge.isRootEdge())
            return logP;
        

        if (edge.parentNode.isReassortment()) {

            BitSet segsToAddLeft = (BitSet) edge.parentNode.getParentEdges().get(0).hasSegments.clone();
            BitSet segsToAddRight = (BitSet) edge.parentNode.getParentEdges().get(1).hasSegments.clone();
            
            for (int segIdx=segsToAdd.nextSetBit(0); segIdx != -1;
                    segIdx=segsToAdd.nextSetBit(segIdx+1)) {
            	
            	if (segsToAddLeft.cardinality()==0 && segsToAddRight.cardinality()==0) { // randomly pick a path when both are empty
                    if (Randomizer.nextBoolean())
                        segsToAddLeft.set(segIdx);
                    else
                        segsToAddRight.set(segIdx);

                    logP += Math.log(0.5);
            	}else if (segsToAddLeft.get(0) && segsToAddRight.cardinality()==1) { // has to follow the core genome when there is already a plasmid transferred
            		segsToAddLeft.set(segIdx);
            	}else if (segsToAddRight.get(0) && segsToAddLeft.cardinality()==1) { // has to follow the core genome when there is already a plasmid transferred
                    segsToAddRight.set(segIdx);
            	}else if (segsToAddLeft.cardinality()>1 && segsToAddRight.cardinality()==1) { // has to follow the path with more plasmids already there as only one plasmid is transferred at a time
            		segsToAddLeft.set(segIdx);
            	}else if (segsToAddRight.cardinality()>1 && segsToAddLeft.cardinality()==1) { // has to follow the path with more plasmids already there as only one plasmid is transferred at a time
            		segsToAddRight.set(segIdx);
            	}else if (segsToAddLeft.cardinality()==1 && segsToAddRight.cardinality()==1) { // as neither path has the core genome, just randomly pick a path
                    if (Randomizer.nextBoolean())
                        segsToAddLeft.set(segIdx);
                    else
                        segsToAddRight.set(segIdx);
                    logP += Math.log(0.5);
            	}else if (segsToAddLeft.cardinality()==1 && segsToAddRight.cardinality()==0) { // left path has one, right has 0, doesn't matter where core goes
                    if (Randomizer.nextDouble()<prob) {
                    	segsToAddRight.set(segIdx);
                        logP += Math.log(prob);
                    }else {
                    	segsToAddLeft.set(segIdx);
                        logP += Math.log(1-prob);
                    }
            	}else if (segsToAddRight.cardinality()==1 && segsToAddLeft.cardinality()==0) { // right path has one, left has 0, doesn't matter where core goes
                    if (Randomizer.nextDouble()<prob) {
                    	segsToAddLeft.set(segIdx);
                        logP += Math.log(prob);
                    }else {
                    	segsToAddRight.set(segIdx);
                        logP += Math.log(1-prob);
                    }
            	}else if (segsToAddLeft.cardinality()>1 && segsToAddRight.cardinality()==0) { // left path has 2+, right has 0, core goes left, otherwise random
            		if (segIdx==0) {
                    	segsToAddLeft.set(segIdx);
            		}else if (Randomizer.nextDouble()<prob) {
            			segsToAddRight.set(segIdx);
                        logP += Math.log(prob);
                    }else {
                    	segsToAddLeft.set(segIdx);
                        logP += Math.log(1-prob);
                    }
            	} else if (segsToAddRight.cardinality()>1 && segsToAddLeft.cardinality()==0) { // right path has 2+, left has 0, core goes right, otherwise random
            		if (segIdx==0) {
            			segsToAddRight.set(segIdx);
            		}else if (Randomizer.nextDouble()<prob) {
            			segsToAddLeft.set(segIdx);
                        logP += Math.log(prob);
                    }else {
                    	segsToAddRight.set(segIdx);
                        logP += Math.log(1-prob);
                    }
            	} else {
            		System.out.println(network);
            		throw new IllegalArgumentException("scenario unknown, should not happen");
            	}
        	}

            logP += addSegmentsToAncestors(edge.parentNode.getParentEdges().get(0), segsToAddLeft);
            logP += addSegmentsToAncestors(edge.parentNode.getParentEdges().get(1), segsToAddRight);
        } else {
            logP += addSegmentsToAncestors(edge.parentNode.getParentEdges().get(0), segsToAdd);
        }

        return logP;
    }
    
    protected BitSet getRandomPlasmid(BitSet sourceSegments) {
        BitSet destSegments = new BitSet();

        destSegments.clear();
        
        int index = 0;
        
        do {
        	index = Randomizer.nextInt(network.getSegmentCount()-1)+1;
        }while (!sourceSegments.get(index));
        
        destSegments.set(index);

        return destSegments;
    }

    protected double getLogUnconditionedPlasmidProb(BitSet sourceSegments) {
        
        int nrPlasmids = sourceSegments.cardinality();
        nrPlasmids -= sourceSegments.get(0) ? 1 : 0;
        
        return Math.log(1.0/nrPlasmids);
    }
    
    


}
