package plasmids.distribution;



import beast.core.CalculationNode;
import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import coalre.distribution.NetworkEvent;
import coalre.network.Network;
import coalre.network.NetworkEdge;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

/**
 * @author Tim Vaughan and Nicola Felix Mueller
 */
public class PlasmidNetworkIntervals extends CalculationNode {
    public Input<Network> networkInput = new Input<>("network",
            "network for which to calculate the intervals", Validate.REQUIRED);
	
	public Input<Function> plasmidTransferRateInput = new Input<>(
	        "plasmidTransferRate",
            "Rate of jumping of an individual plasmid");
	
	

    private Network network;

    private List<NetworkEvent> networkEventList, storedNetworkEventList;

    public boolean eventListDirty = true;
    Function plasmidTransferRate;
    
    boolean hasMultipleRates = false;

    @Override
    public void initAndValidate() {
        network = networkInput.get();

        storedNetworkEventList = new ArrayList<>();
                
        plasmidTransferRate = plasmidTransferRateInput.get();
        
        if (plasmidTransferRate.getDimension()>1)
        	hasMultipleRates = true;
        update();
    }

    public List<NetworkEvent> getNetworkEventList() {
        update();

        return networkEventList;
    }

    void update() {
        if (!eventListDirty)
            return;
        
        networkEventList = network.getNodes().stream().map(n -> {
            NetworkEvent event = new NetworkEvent();
            event.time = n.getHeight();
            event.node = n;
            switch(n.getChildCount()) {
                case 0:
                    event.type = NetworkEvent.NetworkEventType.SAMPLE;
                    break;

                case 1:
                    event.type = NetworkEvent.NetworkEventType.REASSORTMENT;
                    break;

                case 2:
                    event.type = NetworkEvent.NetworkEventType.COALESCENCE;
                    break;

                default:
                    throw new RuntimeException("Network node has illegal number of children.");
            }
            return event;
        }).sorted(Comparator.comparingDouble(e -> e.time)).collect(Collectors.toList());

        int lineages = 0;
        double totalTransferObsProb = 0;

        for (NetworkEvent event : networkEventList) {
            switch(event.type) {
                case SAMPLE:
                    lineages += 1;
                    totalTransferObsProb += getObsProb(event.node.getParentEdges().get(0));
                    break;

                case REASSORTMENT:
                    lineages += 1;
                    
                    totalTransferObsProb -= getObsProb(event.node.getChildEdges().get(0));
                    totalTransferObsProb += getObsProb(event.node.getParentEdges().get(0));
                    totalTransferObsProb += getObsProb(event.node.getParentEdges().get(1));

                    event.segsLeft = (BitSet) event.node.getParentEdges().get(0).hasSegments.clone();
                    event.segsRight = (BitSet) event.node.getParentEdges().get(1).hasSegments.clone();
                    
//                    if (event.segsLeft.cardinality()==0 || event.segsRight.cardinality()==0) {
//                    	System.out.println(event.node.getChildEdges().get(0).childNode.getTaxonLabel());
//
//                    	System.out.println(event.node.getChildEdges().get(0).childNode.getHeight());
//                    	System.out.println(event.node.getChildEdges().get(0).parentNode.getHeight());
//                    	System.out.println(networkInput.get());
//                    	System.exit(0);
//                    }
                    break;

                case COALESCENCE:
                    lineages -= 1;
                    
                    totalTransferObsProb -= getObsProb(event.node.getChildEdges().get(0));
                    totalTransferObsProb -= getObsProb(event.node.getChildEdges().get(1));
                    totalTransferObsProb += getObsProb(event.node.getParentEdges().get(0));
                    break;
            }

            event.lineages = lineages;
            event.totalReassortmentObsProb = totalTransferObsProb;
        }

        eventListDirty = false;
    }
    
    private double getObsProb(NetworkEdge edge) {
    	if (hasMultipleRates) {
    		double prob = 0;
	    	for (int i=1; i < edge.hasSegments.length(); i++) {
	    		if (edge.hasSegments.get(i)) {
	    			prob += plasmidTransferRate.getArrayValue(i-1);
	    		}    			
	    	}
	    	return prob;
    	}
    	
    	if (edge.hasSegments.cardinality()==1)
    		return 0;
    	
        int nrplasmids = edge.hasSegments.cardinality();
        nrplasmids -= edge.hasSegments.get(0) ? 1 : 0;

    	return nrplasmids*plasmidTransferRate.getArrayValue();
    	
    }

    @Override
    protected boolean requiresRecalculation() {
        eventListDirty = true;

        return true;
    }

    @Override
    protected void restore() {
        List<NetworkEvent> tmp = networkEventList;
        networkEventList = storedNetworkEventList;
        storedNetworkEventList = tmp;

        super.restore();
    }

    @Override
    protected void store() {
        storedNetworkEventList.clear();
        storedNetworkEventList.addAll(networkEventList);

        super.store();
    }
}