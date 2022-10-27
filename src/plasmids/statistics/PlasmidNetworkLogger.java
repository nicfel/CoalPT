package plasmids.statistics;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.Tree;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;
import java.util.stream.Collectors;

public class PlasmidNetworkLogger extends BEASTObject implements Loggable {


    public Input<Network> networkInput = new Input<>("network",
            "Network for which to log statistics.",
            Input.Validate.REQUIRED);
   


    Network network;


	@Override
	public void initAndValidate() {
		network = networkInput.get();
	}


    public PlasmidNetworkLogger() { }


    @Override
    public void init(PrintStream out) {
        out.println("#nexus");
        out.println("begin trees;");
    }

    @Override
    public void close(PrintStream out) {
        out.println("End;");
    }

    @Override
    public void log(long sample, PrintStream out) {
        out.println("tree STATE_" + sample + " = " + network.getExtendedNewick(0));
    }


}
