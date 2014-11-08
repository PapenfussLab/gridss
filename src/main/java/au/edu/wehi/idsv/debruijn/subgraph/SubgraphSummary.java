package au.edu.wehi.idsv.debruijn.subgraph;

import java.util.ArrayList;
import java.util.List;

import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;

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
		List<SubgraphSummary> children = new ArrayList<>();
		SubgraphSummary root = start;
		while (root.parent != null) {
			children.add(root);
			root = root.parent;
		}
		for (SubgraphSummary node : children) {
			node.parent = root;
		}
	}
	public static SubgraphSummary merge(SubgraphSummary graph1, SubgraphSummary graph2) {
		graph1 = graph1.getRoot();
		graph2 = graph2.getRoot();
		if (graph1 == graph2) return graph1;
		graph1.parent = graph2;
		graph2.maxAnchor = Math.max(graph2.maxAnchor, graph1.maxAnchor);
		graph2.minAnchor = Math.min(graph2.minAnchor, graph1.minAnchor);
		return graph2;
	}
	private static void addAnchor(SubgraphSummary g, int position) {
		SubgraphSummary rootg = g.getRoot();
		rootg.maxAnchor = Math.max(rootg.maxAnchor, position);
		rootg.minAnchor = Math.min(rootg.minAnchor, position);
	}
	public void addAnchor(int position) {
		addAnchor(this, position);
	}
	public void addAnchor(Integer position) {
		if (position != null) {
			addAnchor(this, position);
		}
	}
	public void addNode(DeBruijnSubgraphNode node) {
		addAnchor(node.getMinReferencePosition());
		addAnchor(node.getMaxReferencePosition());
		addAnchor(node.getMinMatePosition());
		addAnchor(node.getMaxMatePosition());
		// TODO: FIXME: sc with anchor < kmer length will be completely unanchored and subgraph min/max bounds will not be set! 
		// assertion fails
		//assert(node.getMinReferencePosition() != null ||
		//		node.getMaxReferencePosition() != null ||
		//		node.getMinMatePosition() != null || 
		//		node.getMaxMatePosition() != null);
	}
	@Override
	public String toString() {
		return String.format("Subgraph [%d, %d] containing %s", minAnchor, maxAnchor, KmerEncodingHelper.toApproximateString(kmer));
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
