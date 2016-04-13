package au.edu.wehi.idsv.debruijn.positional;

import java.io.File;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.Set;

import au.edu.wehi.idsv.visualisation.PositionalDeBruijnGraphTracker.MemoizationStats;

public abstract class ContigCaller {
	/**
	 * Since reference kmers are not scored, calculating 
	 * highest weighted results in a preference for paths
	 * ending at a RP with sequencing errors over a path
	 * anchored to the reference. 
	 * 
	 * To ensure that the anchored paths are scored higher
	 * than the unanchored paths, paths anchored to the
	 * reference are given a score adjustment larger than
	 * the largest expected score.
	 */
	protected static final int ANCHORED_SCORE = Integer.MAX_VALUE >> 2;
	
	public abstract ArrayDeque<KmerPathSubnode> bestContig(int unprocessedPosition);
	/**
	 * Exports the internal state for debugging purposes
	 * @param file
	 * @throws IOException 
	 */
	public abstract void exportState(File file) throws IOException;
	/**
	 * Called when a node is added to the loaded graph
	 * @param node
	 */
	public void add(KmerPathNode node) { }
	/**
	 * Called when a node is removed from the loaded graph
	 * @param node
	 */
	public void remove(KmerPathNode node) { }
	public void remove(Set<KmerPathNode> keySet) { }

	protected final int maxEvidenceWidth;
	public ContigCaller(int maxEvidenceWidth) {
		this.maxEvidenceWidth = maxEvidenceWidth;
	}
	public abstract boolean sanityCheck();
	public abstract int tracking_memoizedNodeCount();
	public abstract int tracking_frontierSize();
	public abstract MemoizationStats tracking_lastRemoval();
}