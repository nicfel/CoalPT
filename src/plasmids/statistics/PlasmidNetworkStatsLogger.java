package plasmids.statistics;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.stream.Collectors;

public class PlasmidNetworkStatsLogger extends BEASTObject implements Loggable {


    public Input<Network> networkInput = new Input<>("network",
            "Network for which to log statistics.",
            Input.Validate.REQUIRED);
    
    public Input<List<Tree>> segmentTreesInput = new Input<>("segmentTree",
            "Segment tree associated with network.",
            new ArrayList<>());
    
    public Input<RealParameter> plasmidTransferRateInput = new Input<>("plasmidTransferRate",
            "plasmidTransferRate.",
            Input.Validate.OPTIONAL);



    Network network;
    boolean logObservable = false;

    public PlasmidNetworkStatsLogger() { }

    @Override
    public void initAndValidate() {
        network = networkInput.get();
        if(segmentTreesInput.get().size()>0)
        	logObservable = true;
    }

    @Override
    public void init(PrintStream out) {

        String prefix = network.getID() == null ? "networkStat." : network.getID() + ".";

        if (logObservable) {
	        out.print(prefix + "obsHeight\t" +
	                prefix + "obsTotalLength\t");
	        
	        for (int i = 1; i < segmentTreesInput.get().size(); i++) {
	        	out.print("obsPlasmidTransferEvents." + i + "\t");
	        }
	        out.print(prefix + "obsPlasmidTransferEvents.total\t");
         		

        }else {
        	out.print(prefix + "height\t" +
	                prefix + "totalLength\t" +
	                prefix + "PlasmidTransferEvents\t");
        }

    }

    @Override
    public void log(long sample, PrintStream out) {
        if (logObservable){
    		double[] rootHeights = new double[segmentTreesInput.get().size()];
			for (int i = 0; i < segmentTreesInput.get().size(); i++)
				rootHeights[i] = segmentTreesInput.get().get(i).getRoot().getHeight();
			
        	out.print(getTotalHeight(network, rootHeights) + "\t" +
	                getTotalEdgeLength(network, rootHeights) + "\t"); 
	        int totEvents = 0;        
        	
	        for (int i = 1; i < segmentTreesInput.get().size(); i++) {
	        	int plasmidEvents = getReassortmentCount(network, rootHeights, i, plasmidTransferRateInput.get());
	        	totEvents+=plasmidEvents;
	        	out.print(plasmidEvents + "\t");
	        }

	        out.print(totEvents + "\t");
            
        	
        }else{   
        	out.print(getTotalHeight(network) + "\t" +
	                getTotalEdgeLength(network) + "\t" +
	                getReassortmentCount(network) + "\t");
        }
    }

    @Override
    public void close(PrintStream out) {

    }

    public static int getReassortmentCount(Network network) {
        return (int)network.getNodes().stream().filter(NetworkNode::isReassortment).count();
    }

    public static double getTotalEdgeLength(Network network) {
        return network.getEdges().stream().filter(e -> !e.isRootEdge()).
                map(NetworkEdge::getLength).reduce((l1, l2) -> l1+l2).get();
    }

    public static double getTotalHeight(Network network) {
        return network.getRootEdge().childNode.getHeight();
    }
    
    public static int getReassortmentCount(Network network, double[] rootHeights) {
    	double maxHeight = 0.0;
    	for (int i = 0; i < rootHeights.length; i++)
    		if (rootHeights[i] > maxHeight)
    			maxHeight = rootHeights[i];
    	
    	final double finalMaxHeight = maxHeight;

        return (int)network.getNodes().stream()
        		.filter(NetworkNode::isReassortment)
        		.filter(n -> n.getHeight() < finalMaxHeight)
        		.count();
    }
    
    public static int getReassortmentCount(Network network, double[] rootHeights, int segment, RealParameter plasmidTransferRate) {    	
    	final double finalMaxHeight = rootHeights[segment];
    	
    	int count = 0;
    	for (NetworkNode n : network.getNodes()) {
    		if (n.isReassortment() && n.getHeight()<finalMaxHeight){
    			// known main edge
    			if (n.getParentEdges().get(0).hasSegments.get(0) || n.getParentEdges().get(0).hasSegments.cardinality()>1) {
	  				if (n.getParentEdges().get(1).hasSegments.get(segment))
	  					count++;
    			}else if (n.getParentEdges().get(1).hasSegments.get(0) || n.getParentEdges().get(1).hasSegments.cardinality()>1){
	  				if (n.getParentEdges().get(0).hasSegments.get(segment))
	  					count++;
    			}else if (n.getParentEdges().get(0).hasSegments.get(segment) || n.getParentEdges().get(1).hasSegments.get(segment)){
    				// unknown main edge sample from plasmid transfer rates
    				if (plasmidTransferRate.getDimension()==1) {
    					if (Randomizer.nextBoolean())
    						count++;
    				}else {
    					int otherIndex = n.getParentEdges().get(0).hasSegments.get(segment) ? 1 : 0;
    					double ratio = plasmidTransferRate.getArrayValue(segment-1)/
    							(plasmidTransferRate.getArrayValue(segment-1) + 
    									plasmidTransferRate.getArrayValue(n.getParentEdges().get(otherIndex).hasSegments.nextSetBit(0)-1));
    					if (Randomizer.nextDouble()<ratio)
    						count++;
    				}
    			}    			
    		}   		
    		
    	}

        return count;
    }


    public static double getTotalEdgeLength(Network network, double[] rootHeights) {
    	double totalLength = 0.0;
    	List<NetworkEdge> networkEdges = new ArrayList<>(network.getEdges());
        List<NetworkEdge> nonRootEdges = networkEdges.stream()
        		.filter(e -> !e.isRootEdge())
                .collect(Collectors.toList());
        // check for each edge if it has at least one segment for which the root hasn't been reched yet
        for (int i = 0; i < nonRootEdges.size(); i++){
        	double childHeight = nonRootEdges.get(i).childNode.getHeight();
    		final BitSet hasSegment = nonRootEdges.get(i).hasSegments;
    		for (int j = 0; j < rootHeights.length; j++){
    			if (hasSegment.get(j) && childHeight < rootHeights[j]){
    				totalLength += nonRootEdges.get(i).parentNode.getHeight() - childHeight;
    				break;
    			}
    		}
        	
        }        	
    	
        return totalLength;
    }

    public static double getTotalHeight(Network network, double[] rootHeights) {
    	double maxHeight = 0.0;
    	for (int i = 0; i < rootHeights.length; i++)
    		if (rootHeights[i] > maxHeight)
    			maxHeight = rootHeights[i];
        return maxHeight;
    }

}
