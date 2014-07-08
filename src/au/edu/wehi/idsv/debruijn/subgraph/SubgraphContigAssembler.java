package au.edu.wehi.idsv.debruijn.subgraph;

import java.util.LinkedList;
import java.util.List;

public interface SubgraphContigAssembler {
	/**
	 * Extracts all SV-supporting contigs from the subgraph
	 * @return assembled contigs in subgraph
	 */
	List<LinkedList<Long>> assembleContigs();
}
