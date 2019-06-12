package au.edu.wehi.idsv.debruijn.positional;

import au.edu.wehi.idsv.visualisation.PositionalDeBruijnGraphTracker.MemoizationStats;

import java.io.File;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.Set;

public abstract class ContigCaller {
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

	protected final int maxEvidenceSupportIntervalWidth;
	public ContigCaller(int maxEvidenceSupportIntervalWidth) {
		this.maxEvidenceSupportIntervalWidth = maxEvidenceSupportIntervalWidth;
	}
	public abstract boolean sanityCheck();
	public abstract int tracking_memoizedNodeCount();
	public abstract int tracking_frontierSize();
	public abstract MemoizationStats tracking_lastRemoval();
}