package au.edu.wehi.idsv.debruijn;

import java.nio.charset.StandardCharsets;
import java.util.Iterator;
import java.util.List;

import com.google.common.base.Function;
import com.google.common.collect.Iterables;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;

public class PathNode<T extends DeBruijnNodeBase> {
	public int getNodeId() {
		return nodeId;
	}
	public void setNodeId(int nodeId) {
		this.nodeId = nodeId;
	}
	private int nodeId = -1;
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
		// Make our own copies of the positional kmer lists
		for (int i = 0; i < allKmers.size(); i++) {
			allKmers.set(i, Lists.newArrayList(allKmers.get(i)));
		}
		// use already computed stats instead of fully recomputing ourselves
		weight = 0;
		maxKmerWeight = 0;
		for (PathNode<T> pn : path) {
			weight += pn.weight;
			maxKmerWeight = Math.max(maxKmerWeight, pn.maxKmerWeight);
		}
		onKmersChanged(graph);
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
		assert(kmerLength(alternatePath) == length());
		merge(alternatePath, 0, this.length(), 0, graph);
	}
	/**
	 * Merges the kmers on this alternate path into this path
	 * @param alternatePath alternate path to merge
	 * @param alternatePathKmerOffset number of initial kmers on alternate path to skip before merging
	 * @param alternatePathKmers number of kmers on alternate path to merge
	 * @param kmerOffset position to start merging kmers
	 * @param length alternatePath alternate path to merge
	 */
	public void merge(
			Iterable<? extends PathNode<T>> alternatePath,
			int alternatePathKmerOffset,
			int alternatePathKmers,
			int kmerOffset,
			DeBruijnGraphBase<T> graph) {
		assert(kmerLength(alternatePath) >= alternatePathKmers + alternatePathKmerOffset);
		if (alternatePathKmers > length() - kmerOffset) throw new IllegalArgumentException("Alternate path is too long.");
		int i = 0;
		// process all kmers sets on the alternate path
		for (List<Long> nodeKmers :
				Iterables.limit(Iterables.skip(
					Iterables.concat(Iterables.transform(alternatePath, new Function<PathNode<T>, List<List<Long>>>() {
						// flatten to just a sequence of kmer alternatives
						@Override
						public List<List<Long>> apply(PathNode<T> arg) {
							return arg.getPathAllKmers();
						}
					})),
				// grab alternatePathKmers kmers after skipping the first alternatePathKmerOffset kmers
				alternatePathKmerOffset), alternatePathKmers)) { 
			allKmers.get(kmerOffset + i++).addAll(nodeKmers);
		}
		kmersChanged(graph);
	}
	/**
	 * Merges the kmers on this alternate path into this path
	 * @param alternatePath alternate path to merge
	 * @param offset number of initial kmers that the alternate path does not contain.
	 * The first alternate path kmer will be merged into the kmer at offset position.
	 * @param length 
	 */
	public void merge(Iterable<? extends PathNode<T>> alternatePath, int offset, DeBruijnGraphBase<T> graph) {
		if (kmerLength(alternatePath) + offset > path.size()) throw new IllegalArgumentException("Alternate path is too long.");
		int i = 0;
		for (PathNode<T> node : alternatePath) {
			for (List<Long> nodeKmers : node.getPathAllKmers()) {
				allKmers.get(offset + i++).addAll(nodeKmers);
			}
		}
		kmersChanged(graph);
	}
	/**
	 * Iterates over each kmer in the given path
	 * @param it path
	 * @return sequence of primary kmers along the given path
	 */
	public static <T extends DeBruijnNodeBase> Iterator<Long> kmerIterator(Iterable<? extends PathNode<T>> it) {
		return Iterators.concat(Iterators.transform(it.iterator(), new Function<PathNode<T>, Iterator<Long>>() {
				@Override
				public Iterator<Long> apply(PathNode<T> input) {
					return input.path.iterator();
				}
			}));
	}
	/**
	 * Gets the kmer length of the given path
	 * @param it kmer path
	 * @return total kmer length
	 */
	public static <T extends DeBruijnNodeBase> int kmerLength(Iterable<? extends PathNode<T>> it) {
		int len = 0;
		for (PathNode<T> pn : it) {
			len += pn.length();
		}
		return len;
	}
	/**
	 * Gets the total weight of all kmers on the given path
	 * @param path path
	 * @return total weight
	 */
	public static <T extends DeBruijnNodeBase> int kmerTotalWeight(Iterable<? extends PathNode<T>> it) {
		int weight = 0;
		for (PathNode<T> pn : it) {
			weight += pn.getWeight();
		}
		return weight;
	}
	/**
	 * Gets the total weight of all kmers on the given subpath
	 * @param it path
	 * @param kmersToSkip number of kmers
	 * @param kmers number of kmers to calculate weight of
	 * @return total weight
	 */
	public static <T extends DeBruijnNodeBase> int kmerTotalWeight(Iterable<? extends PathNode<T>> it, int kmersToSkip, int kmers, DeBruijnGraphBase<T> graph) {
		int weight = 0;
		for (PathNode<T> pn : it) {
			if (kmersToSkip == 0 && pn.length() <= kmers) {
				// take entire node
				weight += pn.getWeight();
				kmers -= pn.length();
			} else {
				// pro-rata weight according to all the kmers in the positions of interest
				for (long kmer : Iterables.concat(Iterables.limit(Iterables.skip(pn.getPathAllKmers(), kmersToSkip), kmers))) {
					weight += graph.getKmer(kmer).getWeight();
				}
				int kmersSkipped = Math.min(pn.length(), kmersToSkip);
				kmersToSkip -= kmersSkipped;
				int kmersTaken = Math.min(pn.length() - kmersSkipped, kmers);
				kmers -= kmersTaken;
			}
		}
		if (kmers != 0) throw new IllegalArgumentException("Too many kmers requested");
		if (kmersToSkip != 0) throw new IllegalArgumentException("Skipping more kmers than exist on path");
		return weight;
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
	@Override
	public String toString() {
		return toString(null);
	}
	public String toString(DeBruijnGraphBase<T> graph) {
		return String.format("[%4d]%s\t%s", nodeId, printAttributes(), debugPrintPathString(graph));
	}
	protected String printAttributes() {
		int kmerCount = 0;
		for (List<Long> x : allKmers) kmerCount += x.size();
		return String.format(" l=%d\tn=%d\tw=%d", length(), kmerCount - length(), getWeight());
	}
	private String debugPrintPathString(DeBruijnGraphBase<T> graph) {
		if (graph != null) pathKmerString(graph.getK());
		// try to infer k from the states
		int k = 1;
		for (long s : getPath()) {
			while (s >>> (2 * k) > 0) k++;
		}
		return String.format("%s->%s",
				new String(KmerEncodingHelper.encodedToPicardBases(k, getPath().get(0))),
				new String(KmerEncodingHelper.encodedToPicardBases(k, getPath().get(length()-1))));
	}
	public String pathKmerString(int k) {
		return new String(DeBruijnGraphBase.getBaseCalls(this.getPath(), k), StandardCharsets.US_ASCII);
	}
}