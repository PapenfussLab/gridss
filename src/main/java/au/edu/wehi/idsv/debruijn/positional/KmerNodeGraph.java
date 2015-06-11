package au.edu.wehi.idsv.debruijn.positional;

import java.util.Collection;
import java.util.List;

import sun.reflect.generics.reflectiveObjects.NotImplementedException;
import au.edu.wehi.idsv.debruijn.DeBruijnGraph;
import au.edu.wehi.idsv.graph.DirectedAcyclicGraph;

public class KmerNodeGraph implements DeBruijnGraph<KmerNode>, DirectedAcyclicGraph<KmerNode> {
	private int k;
	public KmerNodeGraph(int k) {
		this.k = k;
	}

	@Override
	public int getWeight(KmerNode node) {
		return node.weight();
	}

	@Override
	public List<KmerNode> next(KmerNode node) {
		throw new UnsupportedOperationException();
	}

	@Override
	public List<KmerNode> prev(KmerNode node) {
		throw new UnsupportedOperationException();
	}

	@Override
	public void removeNode(KmerNode node) {
		throw new UnsupportedOperationException();
	}

	@Override
	public void removeEdge(KmerNode source, KmerNode sink) {
		throw new UnsupportedOperationException();
	}

	@Override
	public void addNode(KmerNode node) {
		throw new UnsupportedOperationException();
	}

	@Override
	public void addEdge(KmerNode source, KmerNode sink) {
		throw new UnsupportedOperationException();
	}

	@Override
	public Collection<KmerNode> allNodes() {
		throw new UnsupportedOperationException();
	}

	@Override
	public String toString(Iterable<? extends KmerNode> path) {
		throw new NotImplementedException();
	}

	@Override
	public int getK() {
		return k;
	}

	@Override
	public long getKmer(KmerNode node) {
		return node.kmer();
	}

	@Override
	public boolean isReference(KmerNode node) {
		return node.isReference();
	}

}
