/*
 * Copyright (C) 2022 Nicola Müller <nicola.felix.mueller@gmail.com>
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
 * A rewrite of TreeAnnotator that outputs how often reassortment events happen on trunk branches vs. other branches 
 * @author Nicola Felix Müller <nicola.felix.mueller@gmail.com>
 */
public class PlasmidTransferCount extends ReassortmentAnnotator {

    private enum TrunkDefinition { MostRecentSample, TipDistance }
    
    List<NetworkNode> allTrunkNodes;
    List<Double> leaveDistance;
    List<Boolean> isTrunkNode;

    private static class NetworkAnnotatorOptions {
        File inFile;
        File outFile = new File("transfercount.txt");
        double burninPercentage = 10.0;
        TrunkDefinition trunkDefinition = TrunkDefinition.MostRecentSample;
        double minTipDistance = 2.0;
        int[] removeSegments = new int[0];
        List<File> cladeFiles;


        @Override
        public String toString() {
            return "Active options:\n" +
                    "Input file: " + inFile + "\n" +
                    "Output file: " + outFile + "\n" +
                    "Burn-in percentage: " + burninPercentage + "%\n" +
                    "Definition of the trunk: " + trunkDefinition + "\n" +
            		"minimal distance to a tip to be considered trunk\n" + 
                    "(ignored in MostRecentSample Case): " + minTipDistance;
        }
    }

    public PlasmidTransferCount(NetworkAnnotatorOptions options) throws IOException {

        // Display options:
        System.out.println(options + "\n");
        
        Map<String, Integer> clades = readCladeFiles(options.cladeFiles);
        
        
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
          	ps.print("number\tsegment\tfrom\tto\tfromheight\ttoheight\n");

	        for (Network network : logReader){	
	        	
	        	mapClade(network, clades);
	        	List<NetworkNode> nodes=network.getNodes().stream()
	                    				.filter(e -> e.isReassortment())
	                    				.filter(e -> e.getTypeLabel()!=null)
	                    				.collect(Collectors.toList());
	        	
//	        	System.out.println(network.getExtendedNewick(0));
	        	
	        	for (NetworkNode node : nodes) {
	        		NetworkEdge edge = node.getParentEdges().get(0).hasSegments.get(0) ? node.getParentEdges().get(1) : node.getParentEdges().get(0);
	        		
	        		NetworkNode parent = getCoalParent(edge, edge.hasSegments.nextSetBit(0));
	        		
	        		if (parent!=null)       		
	        			ps.print(counter +"\t" + edge.hasSegments.nextSetBit(0) + 
	        					"\t"+ parent.getTypeLabel() + "\t" + edge.childNode.getTypeLabel() + 
	        					"\t" + parent.getHeight() + "\t" +edge.childNode.getHeight() + "\n");
	        	}
	        	
	        	counter=counter+1;
	        }
	        ps.close();
        }
        System.out.println("\nDone!");
    }
    
    private NetworkNode getCoalParent(NetworkEdge edge, int segment) {
    	if (edge.parentNode==null)
    		return null;
    	
    	if (edge.parentNode.isCoalescence()) {
    		if (edge.parentNode.getChildEdges().get(0).hasSegments.get(segment) && edge.parentNode.getChildEdges().get(1).hasSegments.get(segment))
    			return edge.parentNode;
			else
        		return getCoalParent(edge.parentNode.getParentEdges().get(0), segment);
    	}else{
    		if (edge.parentNode.getParentEdges().get(0).hasSegments.get(segment)) {
        		return getCoalParent(edge.parentNode.getParentEdges().get(0), segment);
    		}else {
    			return getCoalParent(edge.parentNode.getParentEdges().get(1), segment);
    		}
    	}
    }
    
    private void mapClade(Network network, Map<String, Integer> clades) {
    	for (NetworkNode n : network.getLeafNodes()) {    		
    		Integer clade = clades.get(n.getTaxonLabel());
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

        JComboBox<TrunkDefinition> heightMethodCombo = new JComboBox<>(TrunkDefinition.values());

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
                        .addComponent(heightMethodCombo)
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
                        .addComponent(trunkDefinitionLabel)
                        .addComponent(heightMethodCombo,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE))
                .addGroup(layout.createParallelGroup()
                        .addComponent(minTipDistance))
                );

        mainPanel.setBorder(new EtchedBorder());
        cp.add(mainPanel);

        JPanel buttonPanel = new JPanel();

        JButton runButton = new JButton("Analyze");
        runButton.addActionListener((e) -> {
            options.burninPercentage = burninSlider.getValue();
            options.trunkDefinition = (TrunkDefinition)heightMethodCombo.getSelectedItem();
            options.minTipDistance = Double.parseDouble(minTipDistance.getText());
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
                case "-trunkDefinition":
                    if (args.length<=i+1) {
                        printUsageAndError("-trunkDefinition must be either mostRecentSample or minTipDistance.");
                    }

                    try {
                    	if (args[i + 1].equals("mostRecentSample"))
                    		options.trunkDefinition = TrunkDefinition.MostRecentSample;
                    	else if (args[i + 1].equals("minTipDistance"))
                    		options.trunkDefinition = TrunkDefinition.TipDistance;
                    	else
                    		throw new NumberFormatException();

                    } catch (NumberFormatException e) {
                        printUsageAndError("trunkDefinition must be either mostRecentSample or minTipDistance.");
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



                case "-minTipDistance":
                    if (args.length<=i+1) {
                        printUsageAndError("-minTipDistance must be followed by a number.");
                    }

                    try {
                        options.minTipDistance =
                                Double.parseDouble(args[i + 1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("minTipDistance must be a positive number. ");
                     }

                    i += 1;
                    break;
                    
                case "-removeSegments":
                    if (args.length<=i+1) {
                        printUsageAndError("-removeSegments must be followed by at least one number.");
                    }

                    try {
                    	String[] argarray = args[i + 1].split(",");
                    	options.removeSegments = new int[argarray.length];
                    	for (int j = 0; j < argarray.length; j++)
                    		options.removeSegments[j] = Integer.parseInt(argarray[j]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("removeSegments must be an array of integers separated by commas if more than one");
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
            new PlasmidTransferCount(options);

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