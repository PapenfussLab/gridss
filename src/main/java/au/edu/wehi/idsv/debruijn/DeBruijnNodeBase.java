package au.edu.wehi.idsv.debruijn;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import au.edu.wehi.idsv.DirectedEvidence;

import com.google.common.collect.Sets;

public class DeBruijnNodeBase {
	private int nodeWeight = 0;
	private List<DirectedEvidence> supportList = new ArrayList<>();
	public DeBruijnNodeBase(VariantEvidence evidence, int readKmerOffset, ReadKmer kmer) {
		this(kmer.weight, evidence.getDirectedEvidence());
	}
	/**
	 * Creates a node from the given read with the given level of support
	 * @param weight support weight
	 * @param read source read
	 */
	public DeBruijnNodeBase(int weight, DirectedEvidence evidence) {
		if (weight <= 0) throw new IllegalArgumentException("weight must be positive");
		this.nodeWeight = weight;
		supportList.add(evidence);
	}
	/**
	 * Merges the given node into this one
	 * @param node
	 */
	public void add(DeBruijnNodeBase node) {
		this.nodeWeight += node.nodeWeight;
		this.supportList.addAll(node.supportList);
	}
	/**
	 * Reduces the weighting of this node due to removal of a supporting read
	 * @return true if this node is now weightless
	 */
	public void remove(DeBruijnNodeBase node) {
		this.nodeWeight -= node.nodeWeight;
		this.supportList.remove(node.supportList);
	}
	/**
	 * returns the weight of this node
	 * @return weight of this node
	 */
	public int getWeight() {
		return nodeWeight;
	}
	/**
	 * Reads supporting this kmer
	 * @return supporting reads
	 */
	public Set<DirectedEvidence> getSupportingEvidence() {
		return Sets.newHashSet(supportList);
	}
	@Override
	public String toString() {
		return String.format("w=%d, #reads=%d", nodeWeight, getSupportingEvidence().size());
	}
}
