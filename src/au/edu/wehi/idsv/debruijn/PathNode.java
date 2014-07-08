package au.edu.wehi.idsv.debruijn;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import com.google.common.base.Function;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;

public class PathNode<T extends DeBruijnNodeBase> {
	private LinkedList<Long> path;
	private List<LinkedList<Long>> children;
	private int weight;
	private int maxKmerWeight;
	public PathNode(LinkedList<Long> path, DeBruijnGraphBase<T> graph) {
		this.path = path;
		this.children = null;
		this.weight = getPathWeight(path, graph);
		recalcMaxKmerWeight(graph);
		recalc(graph);
	}
	protected void recalc(DeBruijnGraphBase<T> graph) { }
	private void recalcMaxKmerWeight(DeBruijnGraphBase<T> graph) {
		maxKmerWeight = 0;
		List<Iterator<Long>> itList = Lists.newArrayList();
		itList.add(path.iterator());
		if (children != null) {
			for (LinkedList<Long> alt : children) {
				itList.add(alt.iterator());
			}
		}
		for (int i = 0; i < length(); i++) {
			int kmerWeight = 0;
			for (Iterator<Long> it : itList) {
				kmerWeight += graph.getKmer(it.next()).getWeight();
			}
			maxKmerWeight = Math.max(maxKmerWeight, kmerWeight);
		}
	}
	private int getPathWeight(LinkedList<Long> path, DeBruijnGraphBase<T> graph) {
		int weight = 0;
		for (long kmer : path) {
			weight += graph.getKmer(kmer).getWeight();
		}
		return weight;
	}
	public LinkedList<Long> getPath() { return path; }
	/**
	 * Gets kmer paths that have been merged into this path
	 * @return
	 */
	public long getFirst() { return path.getFirst(); }
	public long getLast() { return path.getLast(); }
	public int length() { return path.size(); }
	public int getWeight() { return weight; }
	public int getMaxKmerWeight() { return maxKmerWeight; }
	/**
	 * Kmer paths that have been incorporated into this path
	 * @return alternate kmers merged into this path
	 */
	public List<LinkedList<Long>> getMergedPath() {
		if (children == null) return ImmutableList.of();
		return children;
	}
	/**
	 * Merges the kmers on this alternate path into this path
	 * @param alternatePath alternate path to merge
	 */
	public void merge(LinkedList<Long> alternatePath, DeBruijnGraphBase<T> graph) {
		if (alternatePath.size() != path.size()) throw new IllegalArgumentException("Alternate path must be same length as path");
		if (children == null) children = Lists.newArrayList();
		children.add(alternatePath);
		weight += getPathWeight(alternatePath, graph);
		recalcMaxKmerWeight(graph);
		recalc(graph);
	}
	public static <T extends DeBruijnNodeBase> Iterator<Long> kmerIterator(Iterable<? extends PathNode<T>> it) {
		return Iterators.concat(Iterators.transform(it.iterator(), new Function<PathNode<T>, Iterator<Long>>() {
				@Override
				public Iterator<Long> apply(PathNode<T> input) {
					return input.path.iterator();
				}
			}));
	}
	public boolean contains(long kmer) {
		if (path.contains(kmer)) return true;
		if (children != null) {
			for (LinkedList<Long> alt : children) {
				if (alt.contains(kmer)) return true;
			}
		}
		return false;
	}
}