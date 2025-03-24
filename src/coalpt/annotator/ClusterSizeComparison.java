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
public class ClusterSizeComparison extends ReassortmentAnnotator {
    List<NetworkNode> allTrunkNodes;
    List<Double> leaveDistance;
    List<Boolean> isTrunkNode;
    int stateCount;

    private static class NetworkAnnotatorOptions {
        File inFile;
        File outFile = new File("sizecomp.tsv");
        double burninPercentage = 10.0;
        List<File> cladeFiles;
      

        int segment = 0;

        @Override
        public String toString() {
            return "Active options:\n" +
                    "Input file: " + inFile + "\n" +
                    "Output file: " + outFile + "\n" +
                    "Burn-in percentage: " + burninPercentage + "%\n" +
                    "segment/chromsome or plasmid to use as base tree " + segment + "\n"+
                    "Clade files: " + cladeFiles + "\n";
            
        }
    }

    public ClusterSizeComparison(NetworkAnnotatorOptions options) throws IOException {

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
	      
    	int segmentCount=-1;

        int counter=1;
        
        // compute the pairwise reassortment distances 
        try (PrintStream ps = new PrintStream(options.outFile)) {
        	ps.print("iteration\ttime\tsegment");
        	
        	if (clades.size()>0)
        		ps.print("\tclade");
        	
        	ps.print("\tclusterSize1\tclusterSize2\n");
          	
	        for (Network network : logReader){	  
	        	
	        	segmentCount = network.getSegmentCount();
	        	
	        	if (clades.size()>0)
	        	    mapClade(network, clades);
	        	
	        	Tree segmentTree = new Tree(getTree(network.getRootEdge(), options.segment, Double.POSITIVE_INFINITY, network.getSegmentCount(), clades.size()>0));	        	
	        	// print the header for the file
	        	// calculate the cluster size for each node that has two children with different segments
	        	calculateClusterSize(segmentTree.getRoot(), ps, segmentCount, counter);

	        	counter=counter+1;
	        }
	        ps.close();
        }
        
        System.out.println("\nDone!");
    }
    

	
	private int calculateClusterSize(Node n, PrintStream ps, int nSegments, int counter) {
		int clusterSize = 0;
		if (n.isLeaf()) {
			clusterSize++;
		}else {
			for (Node child : n.getChildren()) {
				clusterSize += calculateClusterSize(child, ps, nSegments, counter);
			}
		}
		n.setMetaData("clusterSize", clusterSize);
		n.metaDataString = n.metaDataString + ",clusterSize=" + n.getMetaData("clusterSize");
		
		// check if this node has two children with different number of segments, if so, print the segment and the 
		// cluster size for each child as well as the node time
		if (n.getChildCount() == 2) {
			for (int i = 1; i < nSegments; i++) {
				if (n.getChild(0).getMetaData("seg"+i)!=n.getChild(1).getMetaData("seg"+i)) {
				 	ps.print(counter + "\t" + n.getHeight() + "\t" + i + "\t");
				 	
				 	
				 	if (n.getChild(0).getMetaData("state")!=null) 
						ps.print(n.getChild(0).getMetaData("state") + "\t");
					
					if ((double) n.getChild(0).getMetaData("seg" + i) == 1.0) {
						ps.print(n.getChild(0).getMetaData("clusterSize") + "\t"
								+ n.getChild(1).getMetaData("clusterSize") + "\n");
					} else {
						ps.print(n.getChild(1).getMetaData("clusterSize") + "\t"
								+ n.getChild(0).getMetaData("clusterSize") + "\n");
					}
				}
			}
		}
		
		return clusterSize;
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
	
    private String getTree(NetworkEdge currentEdge, int segment, double lastCoal, int segmentCount, boolean hasClades) {
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
    	
    	                result.append(getTree(childEdge, segment, currentEdge.childNode.getHeight(), segmentCount, hasClades));
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
    	
    	                result.append(getTree(childEdge, segment, lastCoal, segmentCount, hasClades));
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
			if (hasClades) {				
		        if (currentEdge.childNode.getTypeLabel() != null) {
		        		result.append(",state=\"").append(currentEdge.childNode.getTypeLabel() +"\"");	        		
		        }else {
		        	result.append(",state=\"").append(-1 +"\"");
		        }
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
            new ClusterSizeComparison(options);

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