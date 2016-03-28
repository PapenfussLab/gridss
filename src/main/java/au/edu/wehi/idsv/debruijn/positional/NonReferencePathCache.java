package au.edu.wehi.idsv.debruijn.positional;

import java.util.Collection;
import java.util.SortedSet;

/**
 * Calls optimal contigs from a positional de Bruijn graph
 * 
 * Unlike BestNonReferenceContigCaller, previously calculated
 * paths are retained when calculating successive contigs.
 * 
 * @author Daniel Cameron
 *
 */
public class NonReferencePathCache {
	/**
	 * Terminal nodes sorted by score
	 */
	private SortedSet<MemoizedTraversalNode> byScore;
	private Collection<MemoizedTraversalNode> active;
	/**
	 * Partially computed paths that we cannot process any further
	 */
	private Collection<MemoizedTraversalNode> toUpdate;
	private Collection<KmerPathNode> toProcess;
	/**
	 * Forces a recalculation of all memoized paths involving
	 * the given node
	 * @param node node
	 */
	public void invalidate(KmerPathNode node) {
		// only cache non-reference nodes
	}
	public Contig bestContig() {
		return byScore.isEmpty() ? null : byScore.last();
	}
}
