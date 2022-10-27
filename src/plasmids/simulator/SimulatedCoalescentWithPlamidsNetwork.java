package plasmids.simulator;

import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.math.distributions.MRCAPrior;
import beast.util.Randomizer;
import cern.colt.Arrays;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Comparator;
import java.util.List;

public class SimulatedCoalescentWithPlamidsNetwork extends Network {

    public Input<RealParameter> plasmidTransferRateInput = new Input<>("plasmidTransferRate",
            "Rate of jumping of an individual plasmid", Validate.REQUIRED);

    public Input<PopulationFunction> populationFunctionInput = new Input<>("populationModel",
            "Population model to use.", Validate.REQUIRED);
    
    public Input<List<Tree>> treesInput = new Input<>("tree",
            "multiple trees, the first is assumed to be from the core, all else from plasmids.", new ArrayList<>());

    public Input<Integer> nPlasmidsInput = new Input<>("nPlasmids",
            "Number of segments. Used if no segment trees are supplied.");

    public Input<TraitSet> traitSetInput = new Input<>("traitSet",
            "Trait set used to assign leaf ages.");

    public Input<TaxonSet> taxonSetInput = new Input<>("taxonSet",
            "Taxon set used to define leaves");
    
    public Input<List<TaxonSet>> plasmidTaxonSetInput = new Input<>("plasmidTaxonSet",
            "Taxon set used to define plasmid leaves, needed if not all sequences have all plasmids", new ArrayList<>());

    public Input<Boolean> enableSegTreeUpdateInput = new Input<>("enableSegmentTreeUpdate",
            "If false, segment tree objects won't be updated to agree with simulated " +
                    "network. (Default true.)", true);

    public Input<String> fileNameInput = new Input<>("fileName",
            "Name of file to write simulated network to.");
        

    private PopulationFunction populationFunction;
    private RealParameter plasmidTransferRate;
    

    private int nPlasmids;

    public void initAndValidate() {

        if (treesInput.get().isEmpty()) {
        	nPlasmids = nPlasmidsInput.get();
        }else {
        	nPlasmids = treesInput.get().size()-1;
        }       
               
        segmentCount = nPlasmids+1;

        populationFunction = populationFunctionInput.get();
        plasmidTransferRate = plasmidTransferRateInput.get();

        if (nPlasmids==0) {
            throw new IllegalArgumentException("Need at least one plasmid!");
        }

        // Set up sample nodes:

        List<NetworkNode> sampleNodes = new ArrayList<>();
    

        TaxonSet taxonSet = null;
        if (traitSetInput.get() != null)
            taxonSet = traitSetInput.get().taxaInput.get();
        else if (taxonSetInput.get() != null)
            taxonSet = taxonSetInput.get();
        else if (!treesInput.get().isEmpty())
        	taxonSet = treesInput.get().get(0).getTaxonset();
        else
            throw new IllegalArgumentException("Taxon set must be specified " +
                    "using either taxonSet, traitSet or provided by a segmentTree input.");

        TraitSet traitSet = null;
        if (traitSetInput.get() != null)
            traitSet = traitSetInput.get();
        else if (!treesInput.get().isEmpty())
            traitSet = treesInput.get().get(0).getDateTrait();

        for (int taxonIndex=0; taxonIndex<taxonSet.getTaxonCount(); taxonIndex++) {
            String taxonName = taxonSet.getTaxonId(taxonIndex);

            NetworkNode sampleNode = new NetworkNode();
            sampleNode.setTaxonLabel(taxonName);
            sampleNode.setTaxonIndex(taxonIndex);

            if (traitSet != null)
                sampleNode.setHeight(traitSet.getValue(taxonName));
            else if (!treesInput.get().isEmpty())
                sampleNode.setHeight(treesInput.get().get(0).getNode(taxonIndex).getHeight());
            else
                sampleNode.setHeight(0.0);

            sampleNodes.add(sampleNode);
        }        

        // Perform network simulation:
        simulateNetwork(sampleNodes);

        // Update segment trees:
        if (enableSegTreeUpdateInput.get()) {
            for (int segIdx = 0; segIdx < nPlasmids+1; segIdx++) {
                Tree segmentTree = treesInput.get().get(segIdx);
                updateSegmentTree(segmentTree, segIdx);
                segmentTree.setEverythingDirty(false);
            }
        }
        
//        System.out.println(this);
//        for (int segIdx=0; segIdx<nPlasmids+1; segIdx++)
//            System.out.println(treesInput.get().get(segIdx) +";");
//        
//        System.exit(0);


        // Write simulated network to file if requested
        if (fileNameInput.get() != null) {
            try (PrintStream ps = new PrintStream(fileNameInput.get())) {

                ps.println(toString());

            } catch (FileNotFoundException ex) {
                throw new RuntimeException("Error writing to output file '"
                        + fileNameInput.get() + "'.");
            }
        }

        super.initAndValidate();
    }

    /**
     * Simulate network under coalescent with reassortment model.
     * @param sampleNodes network nodes corresponding to samples.
     */
    public void simulateNetwork(List<NetworkNode> sampleNodes) {

        List<NetworkNode> remainingSampleNodes = new ArrayList<>(sampleNodes);
        List<NetworkEdge> extantLineages = new ArrayList<>();

        remainingSampleNodes.sort(Comparator.comparingDouble(NetworkNode::getHeight));

        double currentTime = 0;
        double timeUntilNextSample;
        do {
            // get the timing of the next sampling event
            if (!remainingSampleNodes.isEmpty()) {
                timeUntilNextSample = remainingSampleNodes.get(0).getHeight() - currentTime;
            } else {
                timeUntilNextSample = Double.POSITIVE_INFINITY;
            }

            // get the current propensities
            int k = extantLineages.size();

            double currentTransformedTime = populationFunction.getIntensity(currentTime);
            double transformedTimeToNextCoal = k>=2 ? Randomizer.nextExponential(0.5*k*(k-1)) : Double.POSITIVE_INFINITY;
                        
            double timeToNextCoal = populationFunction.getInverseIntensity(
                    transformedTimeToNextCoal + currentTransformedTime) - currentTime;
            
            
            double totPlasmidRate = 0.0;
            
            
            for (int i=0; i<nPlasmids;i++) {
            	if (plasmidTransferRate.getDimension()==1)
            		totPlasmidRate+=plasmidTransferRate.getArrayValue();
        		else
            		totPlasmidRate+=plasmidTransferRate.getArrayValue(i);
            }
            
            double timeToNextPlasmidTransfer = k>=1 ? Randomizer.nextExponential(k*totPlasmidRate) : Double.POSITIVE_INFINITY;

            // next event time
            double timeUntilNextEvent = Math.min(timeToNextCoal, timeToNextPlasmidTransfer);
            if (timeUntilNextEvent < timeUntilNextSample) {
                currentTime += timeUntilNextEvent;
                if (timeUntilNextEvent == timeToNextCoal)
                    coalesce(currentTime, extantLineages);
                else
                    transfer(currentTime, extantLineages);
            } else {
                currentTime += timeUntilNextSample;
                sample(remainingSampleNodes, extantLineages);
            }

        }
        while (extantLineages.size() > 1 || !remainingSampleNodes.isEmpty());

        setRootEdge(extantLineages.get(0));
    }

    private void sample(List<NetworkNode> remainingSampleNodes, List<NetworkEdge> extantLineages) {
        // sample the network node
        NetworkNode n = remainingSampleNodes.get(0);

        // Create corresponding lineage
        BitSet hasSegs = new BitSet();
        
        if (plasmidTaxonSetInput.get().size()>0) {
        	hasSegs.set(0, true);
            int c = 1;
            for (TaxonSet t : plasmidTaxonSetInput.get()) {
            	List<String> taxanames = t.asStringList();
            	for (String taxa : taxanames) {
                	if(taxa.contentEquals(n.getTaxonLabel())) {
                		hasSegs.set(c, true);            		
                	}
            	}
            	c++;
            }          	
        }else {
            hasSegs.set(0, nPlasmids+1);
        }
        
        NetworkEdge lineage = new NetworkEdge(null, n, hasSegs);
        extantLineages.add(lineage);
        n.addParentEdge(lineage);

        remainingSampleNodes.remove(0);
    }

    private void coalesce(double coalescentTime, List<NetworkEdge> extantLineages) {
        // Sample the pair of lineages that are coalescing:
        NetworkEdge lineage1 = extantLineages.get(Randomizer.nextInt(extantLineages.size()));
        NetworkEdge lineage2;
        do {
            lineage2 = extantLineages.get(Randomizer.nextInt(extantLineages.size()));
        } while (lineage1 == lineage2);

        // Create coalescent node
        NetworkNode coalescentNode = new NetworkNode();
        coalescentNode.setHeight(coalescentTime)
                .addChildEdge(lineage1)
                .addChildEdge(lineage2);
        lineage1.parentNode = coalescentNode;
        lineage2.parentNode = coalescentNode;

        // Merge segment flags:
        BitSet hasSegments = new BitSet();
        hasSegments.or(lineage1.hasSegments);
        hasSegments.or(lineage2.hasSegments);

        // Create new lineage
        NetworkEdge lineage = new NetworkEdge(null, coalescentNode, hasSegments);
        coalescentNode.addParentEdge(lineage);

        extantLineages.remove(lineage1);
        extantLineages.remove(lineage2);
        extantLineages.add(lineage);
    }

    private void transfer(double reassortmentTime, List<NetworkEdge> extantLineages) {
        NetworkEdge lineage = extantLineages.get(Randomizer.nextInt(extantLineages.size()));

        BitSet hasSegs_left = (BitSet) lineage.hasSegments.clone();
        BitSet hasSegs_right = new BitSet();
        
        double[] cumsum = new double[nPlasmids];
		cumsum[0]=plasmidTransferRate.getArrayValue(0);

        for (int i = 1; i < nPlasmids; i++)
        	if (plasmidTransferRate.getDimension()==1)
        		cumsum[i]=cumsum[i-1] + plasmidTransferRate.getArrayValue();
    		else
        		cumsum[i]=cumsum[i-1] + plasmidTransferRate.getArrayValue(i);
        
        double totrate = cumsum[nPlasmids-1];
        
        double randval = Randomizer.nextDouble();
        
        for (int i = 0; i < nPlasmids; i++)
        	if (randval<=cumsum[i]/totrate) {
        		if (lineage.hasSegments.get(i+1)) {
        			hasSegs_right.set(i+1);
        			hasSegs_left.set(i+1, false);;
        		}
        		break;
        	}else {
        		
        	}

        // Stop here if reassortment event is unobservable
        if (hasSegs_left.cardinality() == 0 || hasSegs_right.cardinality() == 0)
            return;

        // Create reassortment node
        NetworkNode node = new NetworkNode();
        node.setHeight(reassortmentTime).addChildEdge(lineage);

        // Create reassortment lineages
        NetworkEdge leftLineage = new NetworkEdge(null, node, hasSegs_left);
        NetworkEdge rightLineage = new NetworkEdge(null, node, hasSegs_right);
        node.addParentEdge(leftLineage);
        node.addParentEdge(rightLineage);

        extantLineages.remove(lineage);
        extantLineages.add(leftLineage);
        extantLineages.add(rightLineage);
    }
}
