package au.edu.wehi.idsv.debruijn;

import java.nio.charset.StandardCharsets;
import java.util.Iterator;
import java.util.List;

import com.google.common.base.Function;
import com.google.common.collect.Iterables;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;

public class PathNode<T extends DeBruijnNodeBase> {
	private List<Long> path;
	private List<List<Long>> allKmers;
	private int weight;
	private int maxKmerWeight;
	/**
	 * Creates a new path kmer form the given kmer sequence
	 * @param path kmer sequence
	 * @param graph parent graph
	 */
	public PathNode(Iterable<Long> path, DeBruijnGraphBase<T> graph) {
		this.path = Lists.newArrayList(path);
		this.allKmers = Lists.newArrayList(Iterables.transform(path, new Function<Long, List<Long>>() {
			public List<Long> apply(Long arg) {
				return Lists.newArrayList(arg);
			}
		}));
		kmersChanged(graph);
	}
	/**
	 * Merges the given path into a single path node
	 * @param path kmer sequence
	 * @param graph parent graph
	 */
	public PathNode(DeBruijnGraphBase<T> graph, Iterable<? extends PathNode<T>> path) {
		this.path = Lists.newArrayList(Iterables.concat(Iterables.transform(path, new Function<PathNode<T>, List<Long>>() {
			public List<Long> apply(PathNode<T> arg) {
				return arg.path;
			}
		})));
		this.allKmers = Lists.newArrayList(Iterables.concat(Iterables.transform(path, new Function<PathNode<T>, List<List<Long>>>() {
			public List<List<Long>> apply(PathNode<T> arg) {
				return arg.allKmers;
			}
		})));
		// Make our own copy of the positional kmer lists
		for (int i = 0; i < allKmers.size(); i++) {
			allKmers.set(i, Lists.newArrayList(allKmers.get(i)));
		}
		kmersChanged(graph);
	}
	/**
	 * Creates a new path kmer which is a subpath of the given path 
	 * @param unsplit kmer sequence
	 * @param startIndex number of starting kmers of unsplit sequence to skip
	 * @param length length of sequence
	 * @param graph parent graph
	 */
	public PathNode(PathNode<T> unsplit, int startIndex, int length, DeBruijnGraphBase<T> graph) {
		this.path = Lists.newArrayListWithCapacity(length);
		this.allKmers = Lists.newArrayListWithCapacity(length);
		for (int i = startIndex; i < startIndex + length; i++) {
			this.path.add(unsplit.path.get(i));
			this.allKmers.add(Lists.newArrayList(unsplit.allKmers.get(i)));
		}
		kmersChanged(graph);
	}
	private void kmersChanged(DeBruijnGraphBase<T> graph) {
		maxKmerWeight = 0;
		weight = 0;
		for (int i = 0; i < length(); i++) {
			int kmerWeight = 0;
			for (long kmer : allKmers.get(i)) {
				kmerWeight += graph.getKmer(kmer).getWeight();
			}
			weight += kmerWeight;
			maxKmerWeight = Math.max(maxKmerWeight, kmerWeight);
		}
		onKmersChanged(graph);
	}
	protected void onKmersChanged(DeBruijnGraphBase<T> graph) { }
	public List<Long> getPath() { return path; }
	public List<List<Long>> getPathAllKmers() { return allKmers; }
	/**
	 * Gets kmer paths that have been merged into this path
	 * @return
	 */
	public long getFirst() { return path.get(0); }
	public long getLast() { return path.get(length() - 1); }
	public int length() { return path.size(); }
	public int getWeight() { return weight; }
	public int getMaxKmerWeight() { return maxKmerWeight; }
	/**
	 * Merges the kmers on this alternate path into this path
	 * @param alternatePath alternate path to merge
	 */
	public void merge(Iterable<? extends PathNode<T>> alternatePath, DeBruijnGraphBase<T> graph) {
		if (Iterators.size(kmerIterator(alternatePath)) != path.size()) throw new IllegalArgumentException("Alternate path must be same length as path");
		int i = 0;
		for (PathNode<T> node : alternatePath) {
			for (List<Long> nodeKmers : node.getPathAllKmers()) {
				allKmers.get(i++).addAll(nodeKmers);
			}
		}
		kmersChanged(graph);
	}
	public static <T extends DeBruijnNodeBase> Iterator<Long> kmerIterator(Iterable<? extends PathNode<T>> it) {
		return Iterators.concat(Iterators.transform(it.iterator(), new Function<PathNode<T>, Iterator<Long>>() {
				@Override
				public Iterator<Long> apply(PathNode<T> input) {
					return input.path.iterator();
				}
			}));
	}
	public int indexOf(long kmer) {
		for (int i = 0; i < length(); i++) {
			if (allKmers.get(i).contains(kmer)) return i;
		}
		return -1;
	}
	public boolean contains(long kmer) {
		return indexOf(kmer) >= 0;
	}
	private static int nodeCount = 0; // Just used for debugging
	public final int nodeId = nodeCount++; // Just used for debugging
	@Override
	public String toString() {
		return toString(null);
	}
	public String toString(DeBruijnGraphBase<T> graph) {
		int kmerCount = 0;
		for (List<Long> x : allKmers) kmerCount += x.size();
		return String.format("[%5d] l=%d\tn=%d\tw=%d\t%s", nodeId, length(), kmerCount - length(), getWeight(), debugPrintPathString(graph));
	}
	private String debugPrintPathString(DeBruijnGraphBase<T> graph) {
		if (graph != null) return new String(graph.getBaseCalls(path), StandardCharsets.US_ASCII);
		// try to infer k from the states
		int k = 1;
		for (long s : getPath()) {
			while (s >>> (2 * k) > 0) k++;
		}
		return String.format("%s->%s",
				new String(KmerEncodingHelper.encodedToPicardBases(k, getPath().get(0))),
				new String(KmerEncodingHelper.encodedToPicardBases(k, getPath().get(length()-1))));
	}
}