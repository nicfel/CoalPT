/*
 * Copyright (C) 2024 Nicola Müller <nicola.felix.mueller@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package coalpt.annotator;

import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;
import coalre.networkannotator.ReassortmentAnnotator;
import coalre.networkannotator.ReassortmentLogReader;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import java.awt.*;
import java.io.*;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Computes the local branching index from Neher et al. 2014. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4227306/
 * for trees mapped from networks. I.e. first, all segments are mapped onto a particular segment and then the LBI is computed for that segment tree
 * @author Nicola Felix Müller <nicola.felix.mueller@gmail.com>
 */
public class SegmentLBI extends ReassortmentAnnotator {
    List<NetworkNode> allTrunkNodes;
    List<Double> leaveDistance;
    List<Boolean> isTrunkNode;
    int stateCount;

    private static class NetworkAnnotatorOptions {
        File inFile;
        File outFile = new File("tree.trees");
        File outTable;
        double burninPercentage = 10.0;
        List<File> cladeFiles;
        boolean calculateAverageDifference = false;
        

        int segment = 0;
        double tau = 10.0;
        double windowSize = 5.0;


        @Override
        public String toString() {
            return "Active options:\n" +
                    "Input file: " + inFile + "\n" +
                    "Output file: " + outFile + "\n" +
                    "Output Table (if not specified, not table is outputted): " + outTable + "\n" +
                    "Burn-in percentage: " + burninPercentage + "%\n" +
                    "segment/chromsome or plasmid to use as base tree " + segment + "\n"+
                    "tau for LBI calculation: " + tau + "\n"+
                    "Clade files: " + cladeFiles + "\n" +
                    "If true, the average differences between nodes with and without segments is computed: " + calculateAverageDifference + "\n";
            
        }
    }

    public SegmentLBI(NetworkAnnotatorOptions options) throws IOException {

        // Display options:
        System.out.println(options + "\n");
        Map<String, Integer> clades = new HashMap<String, Integer>();
        if (options.cladeFiles!=null) {
            System.out.println("read in clade files: " + options.cladeFiles + "\n");
            clades = readCladeFiles(options.cladeFiles);
            // get all the unique states in clades
            stateCount = clades.values().stream().max(Integer::compare).get()+1;
        }

        // Initialise reader
        ReassortmentLogReader logReader = new ReassortmentLogReader(options.inFile,
                options.burninPercentage);

        System.out.println(logReader.getNetworkCount() + " Networks in file.");

        System.out.println("The first " + logReader.getBurnin() +
                 " (" + options.burninPercentage + "%) ACGs will be discarded " +
                "to account for burnin.");

	      System.out.println("\nWriting output to " + options.outFile.getName()
	      + "...");
	      
    	List<Tree> trees = new ArrayList<>();
    	
    	int segmentCount=-1;

        int counter=1;
        List<Double[]> parameters = new ArrayList<>();
        
        
        // compute the pairwise reassortment distances 
        try (PrintStream ps = new PrintStream(options.outFile)) {
          	ps.print("#NEXUS\n");
          	ps.print("Begin trees;\n");

          	
	        for (Network network : logReader){	  
	        	
	        	segmentCount = network.getSegmentCount();
	        	

	        	if (clades.size()>0)
	        	    mapClade(network, clades);
	        	
	        	
	        	Tree segmentTree = new Tree(getTree(network.getRootEdge(), options.segment, Double.POSITIVE_INFINITY, network.getSegmentCount()));	        	
	        	calculateLBI(segmentTree, options.tau, options.windowSize);
	        	calculateClusterSize(segmentTree);

	    	    if (options.outTable != null) {
		        	trees.add(segmentTree);
	    	    }
//	    	    System.exit(0);
	        	
	        	ps.print("tree STATE_" + counter + " = " + segmentTree
                        + ";\n");
	        	counter=counter+1;
	        }
        	ps.print("End;");

	        ps.close();
        }
        
        
		if (trees.size()>0) {
			System.out.println("\nDone with the trees!");
			int i = 0;
			try (PrintStream ps = new PrintStream(options.outTable)) {
				
				if (options.calculateAverageDifference) { 
					// for each segment at each node, compute what the average difference in LBI or CS for lineages with and without that segment is
					ps.print("iteration\tTime");
					
					for (int j = 0; j < segmentCount; j++) {
						if (j == options.segment)
							continue;
	
						ps.print("\tseg" + j +".CS");
						ps.print("\tseg" + j +".LBI");
						if (stateCount > 1) {
							for (int s = 0; s < stateCount; s++) {
								ps.print("\tseg" + j + ".state" + s + ".CS");
								ps.print("\tseg" + j + ".state" + s + ".LBI");
							}
						}						
					}
					ps.print("\n");
					
					for (Tree t : trees) {
						for (Node n : t.getNodesAsArray()) {
							if (n.isLeaf() || n.isRoot())
								continue;

							ps.print(i + "\t" + n.getHeight());
							
							// for this time, compute the average difference in LBI and CS for lineages with and without segment j

							for (int j = 0; j < segmentCount; j++) {	
								if (j == options.segment)
									continue;
								
								// compare CS and LBI for all lineages with segment j to all without segment j
								List<Double>[] CS  = new ArrayList[2];
								CS[0] = new ArrayList<>();
								CS[1] = new ArrayList<>();
								double[] LBI = new double[2];
								
								int[] linCounts = new int[2];
								
								// also keep track of the state
								List<Double>[][] CS_state  = new ArrayList[2][stateCount];
								for (int k = 0; k < 2; k++) {
									for (int s = 0; s < stateCount; s++) {
										CS_state[k][s] = new ArrayList<>();
									}
								}
								double[][] LBI_state = new double[2][stateCount];
								int[][] linCounts_state = new int[2][stateCount];
								
								for (Node otherNode : t.getNodesAsArray()) {
									if (otherNode.isRoot())
										continue;
									
									if (otherNode.getHeight() <= n.getHeight()
											&& otherNode.getParent().getHeight() > n.getHeight()) {
										if (otherNode.getMetaData("seg" + j) != null) {
											// get the state of other node
											String state = (String) otherNode.getMetaData("state");
											
											if ((double) otherNode.getMetaData("seg" + j) == 1.0) {
												CS[0].add(Math.log((int) otherNode.getMetaData("clusterSize")));
												LBI[0] += (double) otherNode.getMetaData("LBI");
												linCounts[0]++;
												
												if (stateCount > 1) {
													CS_state[0][Integer.parseInt(state)].add(Math.log( (int) otherNode
															.getMetaData("clusterSize")));
													LBI_state[0][Integer.parseInt(state)] += (double) otherNode
															.getMetaData("LBI");
													linCounts_state[0][Integer.parseInt(state)]++;
												}
											} else {
												CS[1].add(Math.log((int) otherNode.getMetaData("clusterSize")));
												LBI[1] += (double) otherNode.getMetaData("LBI");
												linCounts[1]++;
												
												if (stateCount > 1) {
													CS_state[1][Integer.parseInt(state)].add(Math.log( (int) otherNode
															.getMetaData("clusterSize")));
													LBI_state[1][Integer.parseInt(state)] += (double) otherNode
															.getMetaData("LBI");
													linCounts_state[1][Integer.parseInt(state)]++;
												}
											}
																						
										}
									}
								}
								// log standardize the CS after flattening CS
								double mean = 0.0;
								for (int k = 0; k < CS.length; k++) {
									for (Double cs : CS[k]) {
										mean += cs;
									}
                                }
								double std = 0.0;
								mean /= CS[0].size() + CS[1].size();
								for (int k = 0; k < CS.length; k++) {
									for (Double cs : CS[k]) {
										std += Math.pow(cs - mean, 2);
									}
								}
								//compute the mean in CS[0] and CS[1] after log standardizing
								double mean0 = 0.0;
								double mean1 = 0.0;
								for (Double cs : CS[0]) {
									mean0 += (cs - mean) / std;
								}
								for (Double cs : CS[1]) {
									mean1 += (cs - mean) / std;
								}
																
								
								ps.print("\t" + mean0);
								ps.print("\t" + (LBI[0] * linCounts[1]/ (LBI[1] * linCounts[0])));
								
								if (stateCount > 1) {
                                    for (int s = 0; s < stateCount; s++) {
                                    	// log standardize the CS after flattening CS
                                    	double mean_state = 0.0;
                                    	for (int k = 0; k < 2; k++) {
											for (Double cs : CS_state[k][s]) {
												mean_state += cs;
											}
                                    	}
                                    	double std_state = 0.0;
                                    	mean_state /= CS_state[0][s].size() + CS_state[1][s].size();
                                    	for (int k = 0; k < 2; k++) {
                                            for (Double cs : CS_state[k][s]) {
                                                std_state += Math.pow(cs - mean_state, 2);
                                            }
                                    	}
                                    	//compute the mean in CS[0] and CS[1] after log standardizing
                                    	double mean0_state = 0.0;
                                    	for (Double cs : CS_state[0][s]) {
                                            mean0_state += (cs - mean_state) / std_state;
                                                                           	}
                                        ps.print("\t" + mean0_state);
                                        ps.print("\t" + (LBI_state[0][s] * linCounts_state[1][s] / (LBI_state[1][s] * linCounts_state[0][s])));
                                    }
                                }
							}
							ps.print("\n");
          				}
						i++;
					}
				}else {
					ps.print("iteration\tTime\tCS\trelCS.coexist\trelLBI.coexist");
	
					for (int j = 0; j < segmentCount; j++) {
						if (j == options.segment)
							continue;
	
						ps.print("\tseg" + j);
					}
					
					if (stateCount > 1) {
						ps.print("\tstate\trelLBI.coexistClade");					
						ps.print("\trelCS.coexistClade");					
					}
					
					
					ps.print("\n");
					for (Tree t : trees) {
						for (Node n : t.getNodesAsArray()) {
							if (n.isLeaf() || n.isRoot())
								continue;
							
							ps.print(i + "\t" + n.getHeight() +"\t" + n.getMetaData("clusterSize") + "\t" + 
									n.getMetaData("relativeClusterSize") + "\t" + n.getMetaData("relLBI.coexist"));
							
							for (int j = 0; j < segmentCount; j++) {
								if (j == options.segment)
									continue;
								ps.print("\t" + n.getMetaData("seg" + j));
							}
							
							if (stateCount > 1) {
								ps.print("\t" + n.getMetaData("state") +
										"\t" + n.getMetaData("relLBI.coexistClade") +
										"\t" + n.getMetaData("relCS.coexistClade"));
							}
							ps.print("\n");								
						}
						i++;					
					}
				}
			}
		}
        System.out.println("\nDone!");
    }
    

//	private Object[] simulateRandomPlasmid(Network network, Double[] doubles, int i) {
//		// make a new plasmid, plasmid i. First, count how many of each plasmid there are
//		int plasmidCount = 0;
//		int nLeafs = network.getLeafNodes().size();
//		int nPlasmids = network.getSegmentCount();
//		for (NetworkNode n : network.getLeafNodes()) {
//			if (n.getParentEdges().get(0).hasSegments.get(i)) {
//				plasmidCount++;
//			}
//		}
//		
//		// choose plasmidCount unique random numbers between 0 and nLeafs;
//		List<Integer> randomIndices = new ArrayList<>();
//		for (int j = 0; j < plasmidCount; j++) {
//			int random = (int) Math.random() * nLeafs;
//			while (randomIndices.contains(random)) {
//				random = (int) Math.random() * nLeafs;
//			}
//			randomIndices.add(random);
//        }
//		
//		int c = 0;
//		for (NetworkNode n : network.getLeafNodes()) {
//			if (randomIndices.contains(c)) {
//				n.getParentEdges().get(0).hasSegments.set(i);
//			} else {
//				n.getParentEdges().get(0).hasSegments.clear(i);
//			}
//			c++;
//		}
//		
//		
//		// make a new plasmid with the same Ne and plasmid transfer rate
//		
//	}
//

	private List<Double[]> readLogFile (File logFile, int plasmidCount){
    	String line;
    	// skip all lines with #, the first line is the header
    	// look for the header "Ne" and "plasmidTransferRate". If there is more than one plasmidTransferRate, then
    	// they will be denotes using "plasmidTransferRate1", "plasmidTransferRate2", etc.
    	List<Double[]> parameters = new ArrayList<>();
    	List<Integer> indices = new ArrayList<>();
    	try {
        	BufferedReader reader = new BufferedReader(new FileReader(logFile));

        	boolean header = true;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#"))
                    continue;
                
                if (header) {
                	// look for the indices of Ne, and plasmidTransferRate
                	String[] tmp = line.split("\t");
					for (int i = 0; i < tmp.length; i++) {
						if (tmp[i].equals("Ne"))
							indices.add(i);
						if (tmp[i].startsWith("plasmidTransferRate"))
							indices.add(i);
													
					}					                	
                	header = false;
                	continue;
                }
                
                // read in the log files
                String[] tmp = line.split("\t");
                Double[] param = new Double[indices.size()];
				for (int i = 0; i < indices.size(); i++) {
					param[i] = Double.parseDouble(tmp[indices.get(i)]);
				}
				parameters.add(param);
            }                
            reader.close();
            
    	} catch (IOException e) {
    		e.printStackTrace();
    	}
    	return parameters;		
	}


	// Calculate the local branching index for the mapped tree
    private void calculateLBI(Tree segmentTree, double tau, double windowSize) {
		// calculate the LBI using the up down passing from Neher et al. 2014
    	// first, calculate the up passing
    	Node root = segmentTree.getRoot();
    	
    	upPassing(root, tau);
    	downPassing(root, tau, 0.0);		
    	normalizeLBI(root);
		normalizeLBIcoexist(segmentTree, tau, windowSize);
	}

	private double upPassing(Node n, double tau) {
		if (n.isRoot()) {
			double sum = 0.0;
			for (Node child : n.getChildren()) {
				sum += upPassing(child, tau);
			}
			n.setMetaData("up", sum);
			return 0.0;
		}else if (n.isLeaf()) {
			double decay = Math.exp((n.getHeight()-n.getParent().getHeight())/tau);
			double up = tau*(1-decay);
//			System.out.println("up " + up + " decay " + decay + " tau " + tau);
			n.setMetaData("up", 0.0);
            return up;			
		}else {
			double sum = 0.0;
			for (Node child : n.getChildren()) {
				sum += upPassing(child, tau);
			}
			double decay = Math.exp((n.getHeight()-n.getParent().getHeight())/tau);
			
			double up = tau*(1-decay) + decay*sum;
			
			n.setMetaData("up", sum);
			return up;
		}	
	}	
	
	private void downPassing(Node n, double tau, double downContribution) {
		if (!n.isLeaf()) {
			double up = (double) n.getMetaData("up");
			n.setMetaData("LBI", downContribution + up);
			
			for (Node child : n.getChildren()) {
//				System.out.println(child.getHeight() - n.getHeight());
				double decay = Math.exp((child.getHeight() - n.getHeight()) / tau);
				double upChild = tau*(1-decay) + decay*(double) child.getMetaData("up");
				double down = tau*(1-decay) + decay*(downContribution + up - upChild);
				
//				System.out.println("down " + down + " up " + up + " upChild " + upChild + " decay " + decay + " tau " + tau);
				
				
				downPassing(child, tau, down);
			}		
			
		}else {
			n.setMetaData("LBI", downContribution);
		}
		return;		
	}

	private void normalizeLBI(Node root) {
		double maxLBI = getMaxLBI(root);
		normalizeLBIabsolute(root, maxLBI);	
	}
	
	private void normalizeLBIabsolute(Node n, double maxLBI) {
		double normLBI = (double) n.getMetaData("LBI")/maxLBI;
		n.metaDataString = n.metaDataString + ",relLBI.abs=" + normLBI;
		n.setMetaData("relLBI.abs", normLBI);
		if (!n.isLeaf()) {
			for (Node child : n.getChildren()) {
				normalizeLBIabsolute(child, maxLBI);
			}
		}
	}
	
	private void normalizeLBIcoexist(Tree t, double tau) {
		for (Node n : t.getNodesAsArray()) {
			List<Double> ListLBI = new ArrayList<Double>();
			if (n.isLeaf())
				continue;
			
			for (Node otherNode : t.getNodesAsArray()) {
                if (otherNode.isRoot())
                    continue;			
                
                if (otherNode.getHeight()<=n.getHeight() && otherNode.getParent().getHeight()>n.getHeight()) {
                	// compute the LBI of the other lineage based on the time of the current node
                	// get the up value for that node
                	double upChild = (double) otherNode.getMetaData("up");
                	// compute the exponential decay from below until the node
                	double decayBelow = Math.exp((otherNode.getHeight()-n.getHeight())/ tau);
                	// compute the exponential decay from above until the node
                	double decayAbove = Math.exp((n.getHeight()-otherNode.getParent().getHeight())/ tau);
                	// compute the total exponential decay from below to above
                	double decay = Math.exp(otherNode.getHeight()-otherNode.getParent().getHeight() / tau);
                	// compute the contribution from below, treating this time as a single child node
                	double up = tau*(1-decayBelow) + decayBelow*upChild;
                	// get the contribution from above which is just the LBI of the parent node
                	double LBIParent = (double) otherNode.getParent().getMetaData("LBI");
                	// discount the total up contribution of this branch from the LBI
                	double downParent = LBIParent - tau*(1-decayAbove) - decay*up;                	
                	// compute the exponential decay of the down message
                	double down = tau*(1-decayAbove) + decayAbove*downParent;                	                	                	
                	// compute the current LBI assuming this is a single child node
                	double current = down + up;
                	                	
					ListLBI.add(current);
                }
            }
			
			if (ListLBI.size()>2) {
				// get the minimum and maximum entry in the ListLBI
				double minLBI = ListLBI.stream().min(Double::compare).get();
				double maxLBI = ListLBI.stream().max(Double::compare).get();
			
				double normLBI = maxLBI==0.0 ? 0.0 : ((double) n.getMetaData("LBI") - minLBI) / (maxLBI - minLBI);
				n.metaDataString = n.metaDataString + ",relLBI.coexist=" + normLBI;
				n.setMetaData("relLBI.coexist", normLBI);
			}else {
				n.metaDataString = n.metaDataString + ",relLBI.coexist=" + 0.0;
				n.setMetaData("relLBI.coexist", 0.0);

			}
		}
		
		if (stateCount>1) {
			// in this case, normalize the LBI by all coexisting lineages, but only in the same clade
			for (Node n : t.getNodesAsArray()) {

				double maxLBI = 0.0;
				String state = ((String) n.getMetaData("state"));
				if (n.isLeaf())
					continue;
				for (Node otherNode : t.getNodesAsArray()) {
					String otherState = ((String) otherNode.getMetaData("state"));
					if (otherNode.isRoot() || !state.equals(otherState))
						continue;
						
					if (otherNode.getHeight() <= n.getHeight() && otherNode.getParent().getHeight() > n.getHeight()) {
						maxLBI = Math.max(maxLBI, (double) otherNode.getMetaData("LBI"));
					}
				}
				double normLBI = maxLBI==0.0 ? 0.0 : (double) n.getMetaData("LBI") / maxLBI;
				n.metaDataString = n.metaDataString + ",relLBI.coexistClade=" + normLBI;
				n.setMetaData("relLBI.coexistClade", normLBI);

			}

		}
		
	}

	
	private void normalizeLBIcoexist(Tree t, double tau, double windowSize) {
		for (Node n : t.getNodesAsArray()) {
			List<Double> ListLBI = new ArrayList<Double>();
			if (n.isLeaf())
				continue;
			
			for (Node otherNode : t.getNodesAsArray()) {
                if (otherNode.isRoot() || n.isLeaf())
                    continue;			
                
                if (otherNode.getHeight()<=(n.getHeight()+windowSize) && otherNode.getHeight()>(n.getHeight()-windowSize)) {
                	ListLBI.add((double) otherNode.getMetaData("LBI"));
                }
            }
			
			if (ListLBI.size()>2) {
				// get the minimum and maximum entry in the ListLBI
				double minLBI = 0.0;//ListLBI.stream().min(Double::compare).get();
				double maxLBI = ListLBI.stream().max(Double::compare).get();
			
				double normLBI = maxLBI==0.0 ? 0.0 : ((double) n.getMetaData("LBI") - minLBI) / (maxLBI - minLBI);
				n.metaDataString = n.metaDataString + ",relLBI.coexist=" + normLBI;
				n.setMetaData("relLBI.coexist", normLBI);
				
				// check if the relLbi.cpexist is larger than the relLBI.abs
				if (normLBI < (double) n.getMetaData("relLBI.abs")) {
					System.out.println(n.getMetaData("relLBI.abs") + " " + n.getMetaData("relLBI.coexist"));	
					System.out.println(ListLBI);
					System.out.println(maxLBI);
				}
				
			}else {
				n.metaDataString = n.metaDataString + ",relLBI.coexist=" + 0.0;
				n.setMetaData("relLBI.coexist", 0.0);

			}
		}
		
		if (stateCount>1) {
			// in this case, normalize the LBI by all coexisting lineages, but only in the same clade
			for (Node n : t.getNodesAsArray()) {

				List<Double> ListLBI = new ArrayList<Double>();
				
				String state = ((String) n.getMetaData("state"));
				if (n.isLeaf())
					continue;
				
				for (Node otherNode : t.getNodesAsArray()) {
					String otherState = ((String) otherNode.getMetaData("state"));
					if (otherNode.isRoot() || !state.equals(otherState) || otherNode.isLeaf())
						continue;
						
					if (otherNode.getHeight() <= (n.getHeight()+windowSize) 
							&& otherNode.getHeight() > (n.getHeight()-windowSize)) {
						ListLBI.add((double) otherNode.getMetaData("LBI"));
					}
				}
				
				if (ListLBI.size()>2) {

					double minLBI = 0.0;//ListLBI.stream().min(Double::compare).get();
					double maxLBI = ListLBI.stream().max(Double::compare).get();
					
					double normLBI = maxLBI==0.0 ? 0.0 : ((double) n.getMetaData("LBI") - minLBI) / (maxLBI - minLBI);

					n.metaDataString = n.metaDataString + ",relLBI.coexistClade=" + normLBI;
					n.setMetaData("relLBI.coexistClade", normLBI);
				}else {
					n.metaDataString = n.metaDataString + ",relLBI.coexistClade=" + 0.0;
					n.setMetaData("relLBI.coexistClade", 0.0);
				}

			}

		}
		
	}



	private double getMaxLBI(Node n) {
		double max = (double) n.getMetaData("LBI");
		if (n.isLeaf()) {
            return (double) n.getMetaData("LBI");
		}else {
			for (Node child : n.getChildren()) {
                double childLBI = getMaxLBI(child);
                if (childLBI>max)
                    max = childLBI;
            }
			return max;
		}
		
	}
	
	
	private void calculateClusterSize(Tree tree) {
		calculateClusterSize(tree.getRoot());
		computeRelativeClusterSize(tree);		
	}
	
	
	private int calculateClusterSize(Node n) {
		int clusterSize = 0;
		if (n.isLeaf()) {
			clusterSize++;
		}else {
			for (Node child : n.getChildren()) {
				clusterSize += calculateClusterSize(child);
			}
		}
		n.setMetaData("clusterSize", clusterSize);
		n.metaDataString = n.metaDataString + ",clusterSize=" + n.getMetaData("clusterSize");

		return clusterSize;
	}
	
	
	private void computeRelativeClusterSize(Tree tree) {
	    // loop over all internal nodes of the tree,
		// then get all co-existing edges to the timing of that node
		// and compute the relative cluster size from the node and all co-existing edges
		Node[] nodes = tree.getNodesAsArray();
		for (Node n : nodes) {
			if (n.isLeaf())
				continue;
			
			int maxClusterSize = 0;
			List<Double> coexistClusterSize = new ArrayList<>();
			double mean = 0.0;

			for (Node otherNode : nodes) {
			    if (otherNode.isRoot())
			        continue;			

			    if (otherNode.getHeight() <= n.getHeight() && otherNode.getParent().getHeight() > n.getHeight()) {
			        int clusterSize = (int) otherNode.getMetaData("clusterSize");
			        maxClusterSize = Math.max(maxClusterSize, clusterSize);
			        double logClusterSize = Math.log(clusterSize);
			        coexistClusterSize.add(logClusterSize);
			        mean += logClusterSize;
			    }
			}

			mean /= coexistClusterSize.size();

			// Compute the standard deviation
			double std = 0.0;
			for (Double cs : coexistClusterSize) {
			    std += Math.pow(cs - mean, 2);
			}
			std = Math.sqrt(std / coexistClusterSize.size());
			
			if (coexistClusterSize.size()>4 && (int) n.getMetaData("clusterSize")>2) {
				double relativeClusterSize = (Math.log((int) n.getMetaData("clusterSize")) - mean)/std; 
				n.setMetaData("relativeClusterSize",  relativeClusterSize);
				n.metaDataString = n.metaDataString + ",relativeClusterSize=" + n.getMetaData("relativeClusterSize");
			}
		}
		
		if (stateCount>1) {
			// in this case, normalize the LBI by all coexisting lineages, but only in the same clade
			for (Node n : tree.getNodesAsArray()) {

				int maxClusterSize = 0;
				List<Double> coexistClusterSize = new ArrayList<>();
				double mean = 0.0;
				
				String state = ((String) n.getMetaData("state"));
				if (n.isLeaf())
					continue;
				for (Node otherNode : tree.getNodesAsArray()) {
					String otherState = ((String) otherNode.getMetaData("state"));
					if (otherNode.isRoot() || !state.equals(otherState))
						continue;
						
					if (otherNode.getHeight() <= n.getHeight() && otherNode.getParent().getHeight() > n.getHeight()) {
				        int clusterSize = (int) otherNode.getMetaData("clusterSize");
				        maxClusterSize = Math.max(maxClusterSize, clusterSize);
				        double logClusterSize = Math.log(clusterSize);
				        coexistClusterSize.add(logClusterSize);
				        mean += logClusterSize;
					}
				}
				
				mean /= coexistClusterSize.size();

				// Compute the standard deviation
				double std = 0.0;
				for (Double cs : coexistClusterSize) {
				    std += Math.pow(cs - mean, 2);
				}
				std = Math.sqrt(std / coexistClusterSize.size());
				
				if (coexistClusterSize.size()>4 && (int) n.getMetaData("clusterSize")>2) {
					double relativeClusterSize = (Math.log((int) n.getMetaData("clusterSize")) - mean)/std; 
					n.setMetaData("relCS.coexistClade",  relativeClusterSize);
					n.metaDataString = n.metaDataString + ",relCS.coexistClade=" + n.getMetaData("relativeClusterSize");
				}
			}

		}
	}
	

 	private void mapClade(Network network, Map<String, Integer> clades) {
    	for (NetworkNode n : network.getLeafNodes()) {
    		Integer clade = clades.get(n.getTaxonLabel());
    		n.setTypeLabel(Integer.toString(clade));
    		mapCladesOnNetwork(n.getParentEdges().get(0), clade);
    	}		
	}
 	
    
    private void mapCladesOnNetwork(NetworkEdge e, Integer clade) {
    	if (e.isRootEdge())
    		return;
    	if (e.parentNode.getTypeLabel()==null) {
        	e.parentNode.setTypeLabel(Integer.toString(clade));
   		
        	for (NetworkEdge enew : e.parentNode.getParentEdges())
        		if (enew.hasSegments.get(0))
        			mapCladesOnNetwork(enew, clade);
    	}
    }

	public HashMap<String, Integer> readCladeFiles(List<File> cladeFiles) throws IOException {   	
    	
    	HashMap<String, Integer> cladeMap = new HashMap<String, Integer>();
    	
    	int c=0;
    	
    	for (File f : cladeFiles) {
    		BufferedReader reader = new BufferedReader(new FileReader(f));
    		String line = reader.readLine();
    		while (line!=null) {    			
    			String[] tmp = line.split("\\s+");
    			if (tmp.length==0)
    				break;
    			cladeMap.put(tmp[0], c);
    			
    			line = reader.readLine();    			
    		}
    		c++;
    	}
    	return cladeMap;
    }
	
    private String getTree(NetworkEdge currentEdge, int segment, double lastCoal, int segmentCount) {
        StringBuilder result = new StringBuilder();

        if (!currentEdge.childNode.isLeaf()) {
        	if (currentEdge.childNode.isCoalescence() && 
        			currentEdge.childNode.getChildEdges().get(0).hasSegments.get(segment) && 
        			currentEdge.childNode.getChildEdges().get(1).hasSegments.get(segment)) {
        		
                result.append("(");

                boolean isFirst = true;
                for (NetworkEdge childEdge : currentEdge.childNode.getChildEdges()) {
                	if (childEdge.hasSegments.get(segment)) {
    	                if (isFirst)
    	                    isFirst = false;
    	                else
    	                    result.append(",");
    	
    	                result.append(getTree(childEdge, segment, currentEdge.childNode.getHeight(), segmentCount));
                	}
                }

                result.append(")");

    		}else {
                boolean isFirst = true;
                for (NetworkEdge childEdge : currentEdge.childNode.getChildEdges()) {
                	if (childEdge.hasSegments.get(segment)) {
    	                if (isFirst)
    	                    isFirst = false;
    	                else
    	                    result.append(",");
    	
    	                result.append(getTree(childEdge, segment, lastCoal, segmentCount));
                	}
                }
    		}

        	
        }
        
        
        if (currentEdge.childNode.isLeaf() || (currentEdge.childNode.isCoalescence() &&
        		currentEdge.childNode.getChildEdges().get(0).hasSegments.get(segment) && currentEdge.childNode.getChildEdges().get(1).hasSegments.get(segment))) {
	        if (currentEdge.childNode.getTaxonLabel() != null)
	            result.append(currentEdge.childNode.getTaxonLabel());
	
	        result.append("[&");
            result.append("segsCarried=").append(currentEdge.hasSegments.cardinality());
            for (int i=0; i<segmentCount; i++) {
            	if (currentEdge.hasSegments.get(i))
            		result.append(",seg").append(i).append("=1");
                else
                	result.append(",seg").append(i).append("=0");
            }
//	        result.append("segments=").append(currentEdge.hasSegments);
	        if (currentEdge.childNode.getTypeLabel() != null) {
	        		result.append(",state=\"").append(currentEdge.childNode.getTypeLabel() +"\"");	        		
	        }else {
	        	result.append(",state=\"").append(-1 +"\"");
	        }

	        result.append("]");
	
	        if (lastCoal != Double.POSITIVE_INFINITY)
	            result.append(":").append(lastCoal - currentEdge.childNode.getHeight());
	        else
	            result.append(":0.0");
	        
        }

        return result.toString();
    }
	
    /**
     * Use a GUI to retrieve ACGAnnotator options.
     *
     * @param options options object to populate using GUI
     * @return true if options successfully collected, false otherwise
     */
    private static boolean getOptionsGUI(NetworkAnnotatorOptions options) {

        boolean[] canceled = {false};

        JDialog dialog = new JDialog((JDialog)null, true);
        dialog.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        dialog.setLocationRelativeTo(null);
        dialog.setTitle("Plasmid Tree Mapper");

        JLabel logFileLabel = new JLabel("Reassortment Network log file:");
        JLabel outFileLabel = new JLabel("Output file:");
        JLabel burninLabel = new JLabel("Burn-in percentage:");
        JLabel chromosomeIndexLabel = new JLabel("Index of the chromosome (or plasmid)\nto output, starts counting at 0:");

        JTextField inFilename = new JTextField(20);
        inFilename.setEditable(false);
        JButton inFileButton = new JButton("Choose File");

        JTextField outFilename = new JTextField(20);
        outFilename.setText(options.outFile.getName());
        outFilename.setEditable(false);
        JButton outFileButton = new JButton("Choose File");

        JTextField chromosomeIndex = new JTextField(20);
        chromosomeIndex.setText(Integer.toString(options.segment));
        chromosomeIndex.setEditable(true);
//        minTipDistance.setEnabled(false);        

        JSlider burninSlider = new JSlider(JSlider.HORIZONTAL,
                0, 100, (int)(options.burninPercentage));
        burninSlider.setMajorTickSpacing(50);
        burninSlider.setMinorTickSpacing(10);
        burninSlider.setPaintTicks(true);
        burninSlider.setPaintLabels(true);
        burninSlider.setSnapToTicks(true);

        Container cp = dialog.getContentPane();
        BoxLayout boxLayout = new BoxLayout(cp, BoxLayout.PAGE_AXIS);
        cp.setLayout(boxLayout);

        JPanel mainPanel = new JPanel();

        GroupLayout layout = new GroupLayout(mainPanel);
        mainPanel.setLayout(layout);
        layout.setAutoCreateGaps(true);
        layout.setAutoCreateContainerGaps(true);

        layout.setHorizontalGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup()
                        .addComponent(logFileLabel)
                        .addComponent(outFileLabel)
                        .addComponent(burninLabel)
                        .addComponent(chromosomeIndexLabel))
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                        .addComponent(inFilename)
                        .addComponent(outFilename)
                        .addComponent(burninSlider)
                        .addComponent(chromosomeIndex))
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                        .addComponent(inFileButton)
                        .addComponent(outFileButton))
                );

        layout.setVerticalGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup()
                        .addComponent(logFileLabel)
                        .addComponent(inFilename,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE)
                        .addComponent(inFileButton))
                .addGroup(layout.createParallelGroup()
                        .addComponent(outFileLabel)
                        .addComponent(outFilename,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE)
                        .addComponent(outFileButton))
                .addGroup(layout.createParallelGroup()
                        .addComponent(burninLabel)
                        .addComponent(burninSlider,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE))
                .addGroup(layout.createParallelGroup()
                        .addComponent(chromosomeIndexLabel)
                        .addComponent(chromosomeIndex))
                );

        mainPanel.setBorder(new EtchedBorder());
        cp.add(mainPanel);

        JPanel buttonPanel = new JPanel();

        JButton runButton = new JButton("Analyze");
        runButton.addActionListener((e) -> {
            options.burninPercentage = burninSlider.getValue();
            options.segment = Integer.parseInt(chromosomeIndex.getText());
            dialog.setVisible(false);
        });
        runButton.setEnabled(false);
        buttonPanel.add(runButton);

        JButton cancelButton = new JButton("Quit");
        cancelButton.addActionListener((e) -> {
            dialog.setVisible(false);
            canceled[0] = true;
        });
        buttonPanel.add(cancelButton);

        JFileChooser inFileChooser = new JFileChooser();
        inFileButton.addActionListener(e -> {
            inFileChooser.setDialogTitle("Select Reassortment Network log file to summarize");
            if (options.inFile == null)
                inFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
            int returnVal = inFileChooser.showOpenDialog(dialog);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                options.inFile = inFileChooser.getSelectedFile();
                inFilename.setText(inFileChooser.getSelectedFile().getName());
                runButton.setEnabled(true);
            }
        });

        JFileChooser outFileChooser = new JFileChooser();
        outFileButton.addActionListener(e -> {
            outFileChooser.setDialogTitle("Select output file name.");
            if (options.inFile != null)
                outFileChooser.setCurrentDirectory(options.inFile);
            else
                outFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));

            outFileChooser.setSelectedFile(options.outFile);
            int returnVal = outFileChooser.showOpenDialog(dialog);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                options.outFile = outFileChooser.getSelectedFile();
                outFilename.setText(outFileChooser.getSelectedFile().getName());
            }
        });

        cp.add(buttonPanel);

        dialog.pack();
        dialog.setResizable(false);
        dialog.setVisible(true);

        return !canceled[0];
    }

    /**
     * Prepare JFrame to which ACGAnnotator output streams will be
     * directed.
     */
    private static void setupGUIOutput() {

        JFrame frame = new JFrame();
        frame.setTitle("Reassortment Event Locator");
        frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);

        JTextArea textArea = new JTextArea(25, 80);
        textArea.setFont(new Font("monospaced", Font.PLAIN, 12));
        textArea.setEditable(false);
        frame.getContentPane().add(new JScrollPane(textArea), BorderLayout.CENTER);

        JButton closeButton = new JButton("Close");
        closeButton.addActionListener(e -> System.exit(0));
        JPanel buttonPanel = new JPanel();
        buttonPanel.add(closeButton);
        frame.getContentPane().add(buttonPanel, BorderLayout.PAGE_END);

        // Redirect streams to output window:
        OutputStream out = new OutputStream() {
            @Override
            public void write(int b) throws IOException {
                SwingUtilities.invokeLater(() -> {
                    if ((char)b == '\r') {
                        int from = textArea.getText().lastIndexOf("\n") + 1;
                        int to = textArea.getText().length();
                        textArea.replaceRange(null, from, to);
                    } else
                        textArea.append(String.valueOf((char) b));
                });
            }
        };

        System.setOut(new PrintStream(out, true));
        System.setErr(new PrintStream(out, true));

        frame.pack();
        frame.setVisible(true);
    }

    public static String helpMessage =
            "PlasmidTreeMapper - Analyzes reassortment networks to map plasmid trees.\n"
                    + "\n"
                    + "Usage: PlasmidTreeMapper [options] <inputLogFile> [outputFile]\n"
                    + "\n"
                    + "Options:\n"
                    + "--------------------------------------------------------------\n"
                    + "-help                    Display this usage information.\n"
                    + "-burnin <percentage>     Specify the percentage of the log file to discard\n"
                    + "                         as burn-in. Default is 10%.\n"
                    + "-chromosomeIndex <index> Specify the index of the chromosome or plasmid to output,\n"
                    + "                         starting from 0. Default is 0.\n"
                    + "-cladeFileInput <file1,file2,...> Input one or more clade files separated by commas.\n"
                    + "-removeSegments <seg1,seg2,...> Specify segments to remove, separated by commas.\n"
                    + "\n"
                    + "If no output file is specified, the default output file name is 'tree.trees'.\n"
                    + "The inputLogFile is mandatory and must be a valid reassortment network log file.";

    /**
     * Print usage info and exit.
     */
    public static void printUsageAndExit() {
        System.out.println(helpMessage);
        System.exit(0);
    }

    /**
     * Display error, print usage and exit with error.
     */
    public static void printUsageAndError(String errMsg) {
        System.err.println(errMsg);
        System.err.println(helpMessage);
        System.exit(1);
    }

    /**
     * Retrieve TrunkReassortment options from command line.
     *
     * @param args command line arguments
     * @param options object to populate with options
     */
    public static void getCLIOptions(String[] args, NetworkAnnotatorOptions options) {
        int i=0;
        while (args[i].startsWith("-")) {
            switch(args[i]) {
                case "-help":
                    printUsageAndExit();
                    break;

                case "-burnin":
                    if (args.length<=i+1)
                        printUsageAndError("-burnin must be followed by a number (percent)");

                    try {
                        options.burninPercentage = Double.parseDouble(args[i+1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("Error parsing burnin percentage.");
                    }

                    if (options.burninPercentage<0 || options.burninPercentage>100) {
                        printUsageAndError("Burnin percentage must be >= 0 and < 100.");
                    }

                    i += 1;
                    break;
                case "-segment":
                    if (args.length<=i+1) {
                        printUsageAndError("-segment must be one of the segments in the network.");
                    }

                    try {
                        options.segment = Integer.parseInt(args[i + 1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("Error parsing burnin percentage.");
                    }

                    i += 1;
                    break;
				case "-tau":
					if (args.length <= i + 1) {
						printUsageAndError("-tau must be a number.");
					}

					try {
						options.tau = Double.parseDouble(args[i + 1]);
					} catch (NumberFormatException e) {
						printUsageAndError("Error parsing tau.");
					}

					i += 1;
					break;
				case "-windowSize":
					if (args.length <= i + 1) {
						printUsageAndError("-windowSize must be a number.");
					}

					try {
						options.windowSize = Double.parseDouble(args[i + 1]);
					} catch (NumberFormatException e) {
						printUsageAndError("Error parsing tau.");
					}

					i += 1;
					break;

				case "-outTable":
					if (args.length <= i + 1) {
						printUsageAndError("-outTable must a file name.");
					}
					
					try {
						options.outTable = new File(args[i + 1]);
					} catch (NumberFormatException e) {
						printUsageAndError("Error parsing outTable.");
					}					

					i += 1;
					break;

                case "-cladeFileInput":
                    if (args.length<=i+1) {
                        printUsageAndError("-cladeFileInput must be one or more filenames that are separeted by commas.");
                    }

                    try {
                    	String[] filename = args[i + 1].split(",");
                    	options.cladeFiles = new ArrayList<>();
                    	for (int j = 0; j < filename.length; j++) {
                    		options.cladeFiles.add(new File(filename[j]));
                    	}

                    } catch (NumberFormatException e) {
                        printUsageAndError("trunkDefinition must be either mostRecentSample or minTipDistance.");
                    }

                    i += 1;
                    break;
                    
                case "-calculateAverageDifference":
					if (args.length <= i + 1) {
						printUsageAndError("-calculateAverageDifference must be either true or false.");
					}
					try {
						options.calculateAverageDifference = Boolean.parseBoolean(args[i + 1]);
					} catch (NumberFormatException e) {
						printUsageAndError("Error parsing calculateAverageDifference.");
					}
					i += 1;
					break;
											

                default:
                    printUsageAndError("Unrecognised command line option '" + args[i] + "'.");
            }

            i += 1;
        }

        if (i >= args.length)
            printUsageAndError("No input file specified.");
        else
            options.inFile = new File(args[i]);

        if (i+1<args.length)
            options.outFile = new File(args[i+1]);
    }

    /**
     * Main method for ACGAnnotator.  Sets up GUI if needed then
     * uses the ACGAnnotator constructor to actually perform the analysis.
     *
     * @param args command line arguments
     */
    public static void main(String[] args) {
    	NetworkAnnotatorOptions options = new NetworkAnnotatorOptions();

        if (args.length == 0) {
            // Retrieve options from GUI:

            try {
                UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
            } catch (ClassNotFoundException | InstantiationException | UnsupportedLookAndFeelException | IllegalAccessException e) {
                Log.warning.println("Error setting cross-platform look and feel.");
            }

            try {
                SwingUtilities.invokeAndWait(() -> {
                    if (!getOptionsGUI(options))
                        System.exit(0);

                    setupGUIOutput();
                });
            } catch (InterruptedException | InvocationTargetException e) {
                e.printStackTrace();
            }


        } else {
            getCLIOptions(args, options);
        }

        // Run ACGAnnotator
        try {
            new SegmentLBI(options);

        } catch (Exception e) {
            if (args.length == 0) {
                JOptionPane.showMessageDialog(null, e.getMessage(),
                        "Error", JOptionPane.ERROR_MESSAGE);
            } else {
                System.err.println("Error: " + e.getMessage());
                e.printStackTrace();
                System.err.println();
                System.err.println(helpMessage);
            }

            System.exit(1);
        }
    }
}