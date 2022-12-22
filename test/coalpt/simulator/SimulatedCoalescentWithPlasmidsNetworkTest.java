package coalpt.simulator;

import coalpt.simulator.SimulatedCoalescentWithPlamidsNetwork;
import coalpt.statistics.PlasmidNetworkStatsLogger;
import coalre.CoalReTestClass;
import coalre.simulator.SimulatedCoalescentNetwork;

import org.junit.Assert;
import org.junit.Test;

import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.coalescent.ConstantPopulation;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import test.beast.beast2vs1.trace.DiscreteStatistics;

import java.util.List;

public class SimulatedCoalescentWithPlasmidsNetworkTest extends CoalReTestClass {

    @Test
    public void testSimulator() {
        Randomizer.setSeed(1);

        TraitSet dateTrait = getContempDateTraitSet(getTaxonSet(10));
        List<Tree> trees = getSegmentTreeObjects(8, dateTrait);

        ConstantPopulation populationFunction = new ConstantPopulation();
        populationFunction.initByName("popSize", new RealParameter("1.0"));

        int N = 10000;
        double[] reassortmentNodeCounts = new double[N];
        double[] networkHeights = new double[N];
        double[] networkLengths = new double[N];

        for (int i = 0; i < N; i++) {

        	SimulatedCoalescentWithPlamidsNetwork network = new SimulatedCoalescentWithPlamidsNetwork();
            network.initByName(
                    "plasmidTransferRate", new RealParameter("0.1"),
                    "populationModel", populationFunction,
                    "traitSet", dateTrait,
                    "tree", trees.get(0),
                    "tree", trees.get(1),
                    "tree", trees.get(2),
                    "tree", trees.get(3),
                    "tree", trees.get(4),
                    "tree", trees.get(5),
                    "tree", trees.get(6),
                    "tree", trees.get(7),
                    "enableSegmentTreeUpdate", false);

            reassortmentNodeCounts[i] = PlasmidNetworkStatsLogger.getReassortmentCount(network);
            networkHeights[i] = PlasmidNetworkStatsLogger.getTotalHeight(network);
            networkLengths[i] = PlasmidNetworkStatsLogger.getTotalEdgeLength(network);
        }

        double meanCount = DiscreteStatistics.mean(reassortmentNodeCounts);
        double meanHeight = DiscreteStatistics.mean(networkHeights);
        double meanLength = DiscreteStatistics.mean(networkLengths);

        System.out.println(meanCount);
        System.out.println(meanHeight);
        System.out.println(meanLength);

        Assert.assertEquals(4.24, meanCount, 0.1);
        Assert.assertEquals(2.42, meanHeight, 0.1);
        Assert.assertEquals(8.11, meanLength, 0.5);
    }
}
