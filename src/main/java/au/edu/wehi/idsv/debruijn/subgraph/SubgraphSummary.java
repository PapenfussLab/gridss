package au.edu.wehi.idsv.debruijn.subgraph;

import java.util.List;

import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;

import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;
import com.google.common.primitives.Ints;

/**
 * De bruijn subgraph
 * 
 * @author Daniel Cameron
 *
 */
public class SubgraphSummary {
	private SubgraphSummary parent = null;
	private int maxAnchor = Integer.MIN_VALUE;
	private int minAnchor = Integer.MAX_VALUE;
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
	public int getMaxAnchor() {
		return maxAnchor;
	}
	/**
	 * Min anchor position for this subgraph
	 * @return
	 */
	public int getMinAnchor() {
		return minAnchor;
	}
	public boolean isAnchored() {
		return minAnchor != Integer.MAX_VALUE;
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
		SubgraphSummary myRoot = getRoot();
		SubgraphSummary graphRoot = graph.getRoot();
		if (myRoot == graphRoot) return false;
		graphRoot.parent = myRoot;
		myRoot.maxAnchor = Math.max(myRoot.maxAnchor, graphRoot.maxAnchor);
		myRoot.minAnchor = Math.min(myRoot.minAnchor, graphRoot.minAnchor);
		myRoot.kmerCount += graphRoot.kmerCount;
		return true;
	}
	private void addAnchor(SubgraphSummary root, Integer position) {
		if (position != null) {
			root.maxAnchor = Math.max(root.maxAnchor, position);
			root.minAnchor = Math.min(root.minAnchor, position);
		}
	}
	public void addNode(DeBruijnSubgraphNode node) {
		SubgraphSummary root = getRoot();
		addAnchor(root, node.getMinReferencePosition());
		addAnchor(root, node.getMaxReferencePosition());
		addAnchor(root, node.getMinMatePosition());
		addAnchor(root, node.getMaxMatePosition());
		root.kmerCount++;
		// TODO: FIXME: sc with anchor < kmer length will be completely unanchored and subgraph min/max bounds will not be set! 
		// assertion fails
		//assert(node.getMinReferencePosition() != null ||
		//		node.getMaxReferencePosition() != null ||
		//		node.getMinMatePosition() != null || 
		//		node.getMaxMatePosition() != null);
	}
	@Override
	public String toString() {
		return String.format("Subgraph [%d, %d] %d kmers including %s", minAnchor, maxAnchor, kmerCount, KmerEncodingHelper.toApproximateString(kmer));
	}
	public static Ordering<SubgraphSummary> ByMinAnchor = new Ordering<SubgraphSummary>() {
		public int compare(SubgraphSummary g1, SubgraphSummary g2) {
			  return Ints.compare(g1.minAnchor, g2.minAnchor);
		  }
	};
	public static Ordering<SubgraphSummary> ByMaxAnchor = new Ordering<SubgraphSummary>() {
		public int compare(SubgraphSummary g1, SubgraphSummary g2) {
			  return Ints.compare(g1.maxAnchor, g2.maxAnchor);
		  }
	};
}
