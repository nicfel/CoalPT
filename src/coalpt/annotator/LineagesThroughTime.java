/*
 * Copyright (C) 2015 Tim Vaughan <tgvaughan@gmail.com>
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
import cern.colt.Arrays;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;
import coalre.networkannotator.ReassortmentAnnotator;
import coalre.networkannotator.ReassortmentLogReader;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import java.awt.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Returns the number of lineages with and without a plasmid over time.
 * @author Nicola Felix Müller <nicola.felix.mueller@gmail.com>
 */
public class LineagesThroughTime extends ReassortmentAnnotator {


    private static class NetworkAnnotatorOptions {
        File inFile;
        File outFile = new File("ltt.txt");
        double burninPercentage = 10.0;
        double[] timePoints;
        List<File> cladeFiles;
        boolean conditionOnChromosome = false;


        @Override
        public String toString() {
            return "Active options:\n" +
                    "Input file: " + inFile + "\n" +
                    "Output file: " + outFile + "\n" +
                    "Time points: " + Arrays.toString(timePoints) + "\n" +
                    "Clade files: " + cladeFiles + "\n" +
                    "Condition on chromosome: " + conditionOnChromosome + "\n" +
                    "Burn-in percentage: " + burninPercentage + "\n";
        }
    }

    public LineagesThroughTime(NetworkAnnotatorOptions options) throws IOException {
        // Display options:
        System.out.println(options + "\n");
        
        Map<String, Integer> clades = new HashMap<String, Integer>();
        if (options.cladeFiles!=null) {
            System.out.println("read in clade files: " + options.cladeFiles + "\n");
            clades = readCladeFiles(options.cladeFiles);
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
	      

        int counter=1;
        // compute the pairwise reassortment distances 
        try (PrintStream ps = new PrintStream(options.outFile)) {
        	boolean first = true;

	        for (Network network : logReader){	
	        	
	        	if (clades.size()>0)
	        	    mapClade(network, clades);
	        	
	        	List<NetworkEdge> edges=network.getEdges().stream()
        								.filter(e -> !e.isRootEdge())
	                    				.collect(Collectors.toList());
	        	
				if (first) {
					ps.print("iteration\ttime");
					first = false;
					for (int i = 0; i < network.getSegmentCount(); i++) {
						ps.print("\tsegmentprop_" + i);
					}
					
					if (clades.size() > 0) {
						for (int i = 0; i < network.getSegmentCount(); i++) {
							for (int j = 0; j < options.cladeFiles.size(); j++) {
								ps.print("\tclade_" + j + "_segmentprop_" + i);
							}

						}
					}
					
					ps.print("\n");
				}
	        	
        		for (int j = 0; j < options.timePoints.length; j++) {
        			double time = options.timePoints[j];
        			
        			// get all edges for which edge.childNode.getHeight() is below time and edge.parentNode.getHeight() is above time
        			double[] segProb = new double[network.getSegmentCount()];    
        		    int nEdges = 0;
        			// loop over all edges
            		for (NetworkEdge edge : edges) {	
            			// check heights
            			if (options.conditionOnChromosome) {
	            			if (edge.childNode.getHeight() < time && edge.parentNode.getHeight() > time && edge.hasSegments.get(0)) {
								// check if edge carries segment i
								for (int k = 0; k < network.getSegmentCount(); k++)
									if (edge.hasSegments.get(k))
										segProb[k]++;
								nEdges++;
							}
            			}else {
	            			if (edge.childNode.getHeight() < time && edge.parentNode.getHeight() > time) {
	            				// check if edge carries segment i
	            				for (int k = 0; k < network.getSegmentCount();k++)
	            					if (edge.hasSegments.get(k))
	            						segProb[k]++;
	            				nEdges++;
	            			}  
            			}
            		}
            		ps.print(counter + "\t" + time);
					for (int k = 0; k < network.getSegmentCount(); k++)
						ps.print("\t" + segProb[k] / nEdges);
					
					
					if (clades.size() > 0) {
						double[][] segCladeProb = new double[network.getSegmentCount()][options.cladeFiles.size()];
						int[] nEdgesClade = new int[options.cladeFiles.size()];
						// loop over all edges
						for (NetworkEdge edge : edges) {
							// check heights
							if (edge.childNode.getHeight() < time && edge.parentNode.getHeight() > time) {
								if (edge.parentNode.getTypeLabel() == null)
									edge.parentNode.setTypeLabel("-1");
								int clade = Integer.parseInt(edge.parentNode.getTypeLabel());	
								if (clade > -1) {
									nEdgesClade[clade]++;
								}										

								// check if edge carries segment i
								for (int k = 0; k < network.getSegmentCount(); k++) {

									if (edge.hasSegments.get(k)) {
										// parse the string from getTypeLabel such that if the Type label is null the clade is -1										
										if (clade > -1) {
											segCladeProb[k][clade]++;
										}

									}
								}
							}
						} 
						for (int k = 0; k < network.getSegmentCount(); k++)
							for (int l = 0; l < options.cladeFiles.size(); l++) {
								ps.print("\t" + segCladeProb[k][l] / nEdgesClade[l]);
							}
					}
					
					ps.print("\n");
        		}  
	        		        	
	        	counter=counter+1;
	        }
	        ps.close();
        }
        System.out.println("\nDone!");
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
        dialog.setTitle("Reassortment Event Trunk Mapper");

        JLabel logFileLabel = new JLabel("Reassortment Network log file:");
        JLabel outFileLabel = new JLabel("Output file:");
        JLabel burninLabel = new JLabel("Burn-in percentage:");
        JLabel trunkDefinitionLabel = new JLabel("Trunk definition:");

        JTextField inFilename = new JTextField(20);
        inFilename.setEditable(false);
        JButton inFileButton = new JButton("Choose File");

        JTextField outFilename = new JTextField(20);
        outFilename.setText(options.outFile.getName());
        outFilename.setEditable(false);
        JButton outFileButton = new JButton("Choose File");

        JTextField minTipDistance = new JTextField(20);
        minTipDistance.setEditable(true);
//        minTipDistance.setEnabled(false);        

        JSlider burninSlider = new JSlider(JSlider.HORIZONTAL,
                0, 100, (int)(options.burninPercentage));
        burninSlider.setMajorTickSpacing(50);
        burninSlider.setMinorTickSpacing(10);
        burninSlider.setPaintTicks(true);
        burninSlider.setPaintLabels(true);
        burninSlider.setSnapToTicks(true);

//        JSlider thresholdSlider = new JSlider(JSlider.HORIZONTAL,
//                0, 100, (int)(options.convSupportThresh));
//        thresholdSlider.setMajorTickSpacing(50);
//        thresholdSlider.setMinorTickSpacing(10);
//        thresholdSlider.setPaintTicks(true);
//        thresholdSlider.setPaintLabels(true);
//        thresholdSlider.setSnapToTicks(true);

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
                        .addComponent(trunkDefinitionLabel))
//                        .addComponent(thresholdLabel)
//                        .addComponent(geneFlowCheckBox))
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                        .addComponent(inFilename)
                        .addComponent(outFilename)
                        .addComponent(burninSlider)
//                        .addComponent(thresholdSlider)
                        .addComponent(minTipDistance))
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
                        .addComponent(trunkDefinitionLabel))
                .addGroup(layout.createParallelGroup()
                        .addComponent(minTipDistance))
                );

        mainPanel.setBorder(new EtchedBorder());
        cp.add(mainPanel);

        JPanel buttonPanel = new JPanel();

        JButton runButton = new JButton("Analyze");
        runButton.addActionListener((e) -> {
            options.burninPercentage = burninSlider.getValue();
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
            "TrunkReassortment - counts how many reassortment events happened on trunk and non-trunk nodes.\n"
                    + "\n"
                    + "Usage: appstore ACGAnnotator [-help | [options] logFile [outputFile]\n"
                    + "\n"
                    + "Option                   Description\n"
                    + "--------------------------------------------------------------\n"
                    + "-help                    Display usage info.\n"
                    + "-trunkDefinition {MostRecentSample, TipDistance} Choose trunk definition method.\n"
                    + "                         (default MostRecentSample)\n"
                    + "-burnin percentage       Choose _percentage_ of log to discard\n"
                    + "                         in order to remove burn-in period.\n"
                    + "                         (Default 10%)\n"
                    + "-minTipDistance     		minimum distance between internal network node\n"
                    + "                         and tip node such that the internal node is considered trunk.\n"
                    + "                         If not  specified, the trunk is any node between samples\n"
                    + "                         height=0 and the root.\n"
                    + "\n"
                    + "If no output file is specified, output is written to a file\n"
                    + "named 'reassortment_distances.txt'.";

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
				case "-conditionOnChromosome":
					if (args.length <= i + 1)
						printUsageAndError("-conditionOnChromosome must be followed by true or false.");
					
					try {
						options.conditionOnChromosome = Boolean.parseBoolean(args[i + 1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("Error parsing conditionOnChromosome.");
                    }
					i += 1;
					break;
					
                case "-timepoints":
					if (args.length <= i + 1)
						printUsageAndError("-timepoints must be followed by 3 numbers start,end,stepsize.");

					String[] timePoints = args[i + 1].split(",");
					// check that there are three values
					if (timePoints.length != 3)
                        printUsageAndError("-timepoints must be followed by 3 numbers start,end,stepsize.");
					
					try {
						double currtime = Double.parseDouble(timePoints[0]);
						double endtime = Double.parseDouble(timePoints[1]);
						double stepsize = Double.parseDouble(timePoints[2]);
						
						int numTimePoints = (int) ((endtime - currtime) / stepsize) + 1;
						options.timePoints = new double[numTimePoints];
						for (int j = 0; j < numTimePoints; j++) {
							options.timePoints[j] = currtime;
							currtime += stepsize;
						}

					} catch (NumberFormatException e) {
						printUsageAndError("Error parsing timepoints.");
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
            new LineagesThroughTime(options);

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