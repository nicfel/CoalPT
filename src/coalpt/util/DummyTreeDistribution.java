package coalpt.util;

import beast.base.core.Input;
import beast.base.evolution.tree.TreeDistribution;
/**
 * Dummy object used only by BEAUti template.
 */
public class DummyTreeDistribution extends TreeDistribution {
    public Input<Boolean> isChromosome = new Input<>("isChromosome", "If yes, this tree is used as the chromosomal tree.", false);
}
