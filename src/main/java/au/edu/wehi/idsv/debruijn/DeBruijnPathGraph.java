package au.edu.wehi.idsv.debruijn;

import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.graph.PathGraph;
import au.edu.wehi.idsv.graph.PathNodeFactory;
import au.edu.wehi.idsv.graph.WeightedSequenceGraphNodeUtil;
import au.edu.wehi.idsv.util.AlgorithmRuntimeSafetyLimitExceededException;
import au.edu.wehi.idsv.visualisation.SubgraphAssemblyAlgorithmTracker;
import au.edu.wehi.idsv.visualisation.TimeoutNodeTraversalTracker;
import com.google.common.base.Predicate;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Queues;

import java.util.Collection;
import java.util.Deque;
import java.util.List;
import java.util.PriorityQueue;

public class DeBruijnPathGraph<T, PN extends DeBruijnPathNode<T>> extends PathGraph<T, PN> implements DeBruijnGraph<PN> {
	public DeBruijnPathGraph(DeBruijnGraph<T> graph, Collection<T> seeds, PathNodeFactory<T, PN> factory, SubgraphAssemblyAlgorithmTracker<T, PN> tracker) {
		super(graph, seeds, factory, tracker);
	}
	public DeBruijnPathGraph(DeBruijnGraph<T> graph, PathNodeFactory<T, PN> factory, SubgraphAssemblyAlgorithmTracker<T, PN> tracker) {
		super(graph, factory, tracker);
	}
	@Override
	public DeBruijnGraph<T> getGraph() {
		return (DeBruijnGraph<T>)super.getGraph();
	}
	@Override
	public boolean isReference(PN node) {
		return node.isReference();
	}
	@Override
	public long getKmer(PN node) {
		throw new UnsupportedOperationException("path nodes contain multiple kmers");
	}
	@Override
	public int getK() {
		return getGraph().getK();
	}
}
