package au.edu.wehi.idsv.debruijn.subgraph;

import java.util.List;

import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;

import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;
import com.google.common.primitives.Longs;

/**
 * De bruijn subgraph
 * 
 * @author Daniel Cameron
 *
 */
public class SubgraphSummary {
	private SubgraphSummary parent = null;
	private long maxPosition = Long.MIN_VALUE;
	private long minPosition = Long.MAX_VALUE;
	private int kmerCount = 0;
	private long kmer;
	public SubgraphSummary(long startingKmer) {
		this.kmer = startingKmer;
	}
	/**
	 * Returns a kmer which is part of this subgraph
	 * @return
	 */
	public long getAnyKmer() {
		return kmer;
	}
	/**
	 * Max anchor position for this subgraph
	 * @return
	 */
	public long getMaxLinearPosition() {
		return maxPosition;
	}
	/**
	 * Max anchor position for this subgraph
	 * @return
	 */
	public long getMinLinearPosition() {
		return minPosition;
	}
	/**
	 * Number of kmers in this subgraph
	 * @return
	 */
	public int getKmerCount() {
		return kmerCount;
	}
	/**
	 * Gets the full subgraph containing this partial subgraph
	 * @return Containing subgraph
	 */
	public SubgraphSummary getRoot() {
		if (parent == null) return this; // we are root
		if (parent.parent != null) collapseToAncestor(this);
		return parent; // we are a child of the root
	}
	/**
	 * Collapses this and all parents so they are all direct children of the root
	 * Amortised cost of this operation is O(1) as all future calls to getRoot()
	 * on all ancestor nodes will immediately return the root
	 */
	private static void collapseToAncestor(SubgraphSummary start) {
		List<SubgraphSummary> children = Lists.newArrayList();
		SubgraphSummary root = start;
		while (root.parent != null) {
			children.add(root);
			root = root.parent;
		}
		for (SubgraphSummary node : children) {
			node.parent = root;
		}
	}
	/**
	 * Merges the given graph into this one
	 * @param graph graph to merge
	 * @return true if the given graph was a distinct graph that merged into this one, false if
	 * both graphs represent the same graph
	 */
	public boolean add(SubgraphSummary graph) {
		SubgraphSummary root = getRoot();
		SubgraphSummary graphRoot = graph.getRoot();
		if (root == graphRoot) return false;
		graphRoot.parent = root;
		root.maxPosition = Math.max(root.maxPosition, graphRoot.maxPosition);
		root.minPosition = Math.max(root.minPosition, graphRoot.minPosition);
		root.kmerCount += graphRoot.kmerCount;
		return true;
	}
	public void addNode(DeBruijnSubgraphNode node) {
		SubgraphSummary root = getRoot();
		root.maxPosition = Math.max(root.maxPosition, node.getMaxLinearPosition());
		root.minPosition = Math.max(root.minPosition, node.getMinLinearPosition());
		root.kmerCount++;
	}
	@Override
	public String toString() {
		return String.format("Subgraph [%d,%d] %d kmers including %s", getMinLinearPosition(), getMaxLinearPosition(), kmerCount, KmerEncodingHelper.toApproximateString(kmer));
	}
	public static Ordering<SubgraphSummary> ByMaxPosition = new Ordering<SubgraphSummary>() {
		public int compare(SubgraphSummary g1, SubgraphSummary g2) {
			  return Longs.compare(g1.maxPosition, g2.maxPosition);
		  }
	};
}
