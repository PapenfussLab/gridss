package au.edu.wehi.idsv.debruijn.positional;

import it.unimi.dsi.fastutil.ints.Int2BooleanSortedMap;
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.NavigableSet;
import java.util.Set;
import java.util.SortedSet;

import au.edu.wehi.idsv.debruijn.DeBruijnGraph;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.graph.DirectedAcyclicGraph;

import com.google.common.base.Function;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeMap;
import com.google.common.collect.TreeRangeSet;

public class KmerAggregrateGraph implements DeBruijnGraph<KmerAggregateNode>, DirectedAcyclicGraph<KmerAggregateNode> {
	private final int k;
	private final Map<Long, RangeMap<Integer, Integer>> kmerIntervalWeights;
	private final Long2ObjectOpenHashMap<Int2BooleanSortedMap> reference = null; // TODO: populate this map
	public KmerAggregrateGraph(int k, Map<Long, SortedSet<KmerSupportNode>> support) {
		this.k = k;
		this.kmerIntervalWeights = transform(support);
		splitIntervals();
	}
	private static Map<Long, RangeMap<Integer, Integer>> transform(Map<Long, SortedSet<KmerSupportNode>> support) {
		Map<Long, RangeMap<Integer, Integer>> result = new HashMap<Long, RangeMap<Integer,Integer>>();
		for (Entry<Long, SortedSet<KmerSupportNode>> entry : support.entrySet()) {
			NavigableSet<KmerAggregateNode> split = KmerAggregateNode.aggregate(entry.getKey(), entry.getValue());
			RangeMap<Integer, Integer> rm = TreeRangeMap.create();
			for (KmerAggregateNode n : split) {
				rm.put(Range.closedOpen(n.getStartPosition(), n.getEndPosition() + 1), n.getWeight());
			}
			result.put(entry.getKey(), rm);
		}
		return result;
	}
	/**
	 * Splits intervals until all previous and next nodes contain full intervals
	 * @param m interval kmer map
	 */
	public void splitIntervals() {
		Set<Long> frontier = new HashSet<Long>(kmerIntervalWeights.keySet());
		while (!frontier.isEmpty()) {
			long kmer = frontier.iterator().next();
			frontier.remove(kmer);
			if (splitIntervals(kmer) > 0) {
				// add adjacent kmers as they may now need to be split
				for (long adjKmer : KmerEncodingHelper.adjacentStates(k, kmer)) {
					frontier.add(adjKmer);	
				}
			}
		}
	}
	private int splitIntervals(long kmer) {
		RangeMap<Integer, Integer> rm = kmerIntervalWeights.get(kmer);
		if (rm == null) return 0;
		Set<Integer> expectedBoundaryValues = new HashSet<Integer>();
		int delta = 1;
		for (long adjKmer : KmerEncodingHelper.nextStates(k, kmer)) {
			RangeMap<Integer, Integer> intervals = kmerIntervalWeights.get(adjKmer);
			if (intervals != null) {
				for (Range<Integer> r : intervals.asMapOfRanges().keySet()) {
					expectedBoundaryValues.add(r.lowerEndpoint() - delta);
					expectedBoundaryValues.add(r.upperEndpoint() - delta);
				}
			}
		}
		delta = -1;
		for (long adjKmer : KmerEncodingHelper.prevStates(k, kmer)) {
			RangeMap<Integer, Integer> intervals = kmerIntervalWeights.get(adjKmer);
			if (intervals != null) {
				for (Range<Integer> r : intervals.asMapOfRanges().keySet()) {
					expectedBoundaryValues.add(r.lowerEndpoint() - delta);
					expectedBoundaryValues.add(r.upperEndpoint() - delta);
				}
			}
		}
		int changes = splitIntervals(rm, expectedBoundaryValues);
		return changes;
	}
	private static int splitIntervals(RangeMap<Integer, Integer> rm, Collection<Integer> expectedBoundaryValues) {
		int changedCount = 0;
		for (int bound : expectedBoundaryValues) {
			Entry<Range<Integer>, Integer> entry = rm.getEntry(bound);
			Range<Integer> currentRange = entry.getKey();
			if (currentRange.lowerEndpoint() != bound) {
				int weight = entry.getValue();
				changedCount++;
				rm.remove(currentRange);
				rm.put(Range.closedOpen(currentRange.lowerEndpoint(), bound), weight);
				rm.put(Range.closedOpen(bound, currentRange.upperEndpoint()), weight);
			}
		}
		return changedCount;
	}
	public List<KmerAggregateNode> next(KmerAggregateNode node) {
		return adj(node, KmerEncodingHelper.nextStates(k, node.getKmer()), 1);
	}
	public List<KmerAggregateNode> prev(KmerAggregateNode node) {
		return adj(node, KmerEncodingHelper.prevStates(k, node.getKmer()), -1);
	}
	private List<KmerAggregateNode> adj(KmerAggregateNode node, long[] states, int delta) {
		List<KmerAggregateNode> adj = new ArrayList<KmerAggregateNode>(4);
		for (long kmer : states) {
			RangeMap<Integer, Integer> intervals = kmerIntervalWeights.get(kmer);
			if (intervals != null) {
				Entry<Range<Integer>, Integer> entry = intervals.getEntry(node.getStartPosition() + delta);
				if (entry != null) {
					// assumes that the intervals have alread been split
					assert(entry.getKey().lowerEndpoint() == node.getStartPosition() + delta);
					assert(entry.getKey().upperEndpoint() == node.getEndPosition() + 1 + delta);
					adj.add(new KmerAggregateNode(kmer, entry.getValue(), entry.getKey().lowerEndpoint(), entry.getKey().upperEndpoint()));
				}
			}
		}
		return adj;
	}
	public List<KmerAggregateNode> nextUnsplit(KmerAggregateNode node) {
		return nextUnsplit(node.getKmer(), node.getStartPosition(), node.getEndPosition());
	}
	public List<KmerAggregateNode> nextUnsplit(long kmer, int start, int end) {
		return adjUnsplit(start, end, KmerEncodingHelper.nextStates(k, kmer), 1);
	}
	public List<KmerAggregateNode> prevUnsplit(KmerAggregateNode node) {
		return prevUnsplit(node.getKmer(), node.getStartPosition(), node.getEndPosition());
	}
	public List<KmerAggregateNode> prevUnsplit(long kmer, int start, int end) {
		return adjUnsplit(start, end, KmerEncodingHelper.prevStates(k, kmer), -1);
	}
	public List<KmerAggregateNode> adjUnsplit(int start, int end, long[] adj, int delta) {
		List<KmerAggregateNode> states = new ArrayList<KmerAggregateNode>(4);
		for (long adjKmer : adj) {
			RangeMap<Integer, Integer> intervals = kmerIntervalWeights.get(adjKmer);
			if (intervals != null) {
				for (Entry<Range<Integer>, Integer> entry : intervals.subRangeMap(Range.closedOpen(start + delta, end + 1 + delta)).asMapOfRanges().entrySet()) {
					states.add(new KmerAggregateNode(adjKmer, entry.getValue(), entry.getKey().lowerEndpoint(), entry.getKey().upperEndpoint()));
				}
			}
		}
		return states;
	}
	/**
	 * Returns all starting nodes for which there is no connected previous kmer position
	 * @return all starting nodes
	 */
	public List<KmerAggregateNode> getStartingNodesUnsplit() {
		List<KmerAggregateNode> startingNodes = new ArrayList<KmerAggregateNode>(4);
		for (Entry<Long, RangeMap<Integer, Integer>> kmerEntry : kmerIntervalWeights.entrySet()) {
			long kmer = kmerEntry.getKey();
			RangeMap<Integer, Integer> intervals = kmerEntry.getValue();
			for (Entry<Range<Integer>, Integer> entry : intervals.asMapOfRanges().entrySet()) {
				int weight = entry.getValue();
				int start = entry.getKey().lowerEndpoint();
				int end = entry.getKey().upperEndpoint() - 1;
				startingNodes.addAll(getStartingNodesUnsplit(kmer, weight, start, end));
			}
		}
		return startingNodes;
	}
	/**
	 * Calculates the subintervals for the given kmer interval for which there is no previous value 
	 * @param kmer kmer
	 * @param weight weight of kmer interval
	 * @param start starting position
	 * @param end ending position
	 * @return kmer nodes containing the disjoint interval subsets for which there is no previous kmer
	 */
	private List<KmerAggregateNode> getStartingNodesUnsplit(long kmer, int weight, int start, int end) {
		// asser the queried value exists as expected
		assert(kmerIntervalWeights.get(kmer) != null);
		assert(kmerIntervalWeights.get(kmer).asMapOfRanges().containsKey(Range.closedOpen(start, end + 1)));
		assert(kmerIntervalWeights.get(kmer).asMapOfRanges().get(Range.closedOpen(start, end + 1)) == weight);
		// find the intervals where there exists a connected prev kmer
		int delta = -1;
		RangeSet<Integer> existing = TreeRangeSet.create();
		for (long adjKmer : KmerEncodingHelper.prevStates(k, kmer)) {
			RangeMap<Integer, Integer> intervals = kmerIntervalWeights.get(adjKmer);
			if (intervals != null) {
				for (Range<Integer> r : intervals.subRangeMap(Range.closedOpen(start + delta, end + 1 + delta)).asMapOfRanges().keySet()) {
					existing.add(Range.closedOpen(r.lowerEndpoint() - delta, r.upperEndpoint() - delta));
				}
			}
		}
		// input [1-9] -> [1-10) interval
		// { [1-2), [4,5), [5, 6), [9,10) } ->
		// { [1-2), [4,6), [9,10) } -> complement()
		// { [2, 4), [6, 9) } ->
		// { [2, 3], [6, 8] } which are the interval for which the input has no prev kmer
		RangeSet<Integer> starting = existing.complement().subRangeSet(Range.closedOpen(start, end + 1));
		List<KmerAggregateNode> startingNodes = new ArrayList<KmerAggregateNode>(starting.asRanges().size());
		for (Range<Integer> r : starting.asRanges()) {
			startingNodes.add(new KmerAggregateNode(kmer, weight, r.lowerEndpoint(), r.upperEndpoint() - 1));
		}
		return startingNodes;
	}
	@Override
	public int getWeight(KmerAggregateNode node) {
		return node.getWeight();
	}
	@Override
	public void removeNode(KmerAggregateNode node) {
		throw new UnsupportedOperationException();
	}
	@Override
	public void removeEdge(KmerAggregateNode source, KmerAggregateNode sink) {
		throw new UnsupportedOperationException();		
	}
	@Override
	public void addNode(KmerAggregateNode node) {
		throw new UnsupportedOperationException();
	}
	@Override
	public void addEdge(KmerAggregateNode source, KmerAggregateNode sink) {
		throw new UnsupportedOperationException();
	}
	@Override
	public Collection<KmerAggregateNode> allNodes() {
		throw new UnsupportedOperationException();
	}
	@Override
	public String toString(Iterable<? extends KmerAggregateNode> path) {
		Iterable<Long> kmerIt = Iterables.transform(path, new Function<KmerAggregateNode, Long>() {
			@Override
			public Long apply(KmerAggregateNode input) {
				return input.getKmer();
			}
		});
		return new String(KmerEncodingHelper.baseCalls(Lists.newArrayList(kmerIt), k), StandardCharsets.US_ASCII);
	}
	@Override
	public boolean isReference(KmerAggregateNode node) {
		Int2BooleanSortedMap lookup = reference.get(node.getKmer());
		// requires sentinel placeholders outside of bounds
		// TODO: no lowerEntry() in fastutils - is heapMap() our only option?
		return lookup.get(lookup.headMap(node.getStartPosition()).firstIntKey());
	}
	@Override
	public long getKmer(KmerAggregateNode node) {
		return node.getKmer();
	}
	@Override
	public int getK() {
		return k;
	}
}
