package au.edu.wehi.idsv.debruijn;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.base.Function;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Ordering;
import com.google.common.collect.Sets;
import com.google.common.primitives.Ints;
import com.google.common.primitives.Longs;

import au.edu.wehi.idsv.visualisation.DeBruijnGraphExporter;
import gnu.trove.map.hash.TLongObjectHashMap;
import gnu.trove.set.TLongSet;
import gnu.trove.set.hash.TLongHashSet;

/**
 * Basic De Bruijn graph implementation
 * @author Daniel Cameron
 *
 * @param <T>
 */
public abstract class DeBruijnGraphBase<T extends DeBruijnNodeBase> implements DeBruijnGraph<T> {
	public static final byte MAX_QUAL_SCORE = 128 - 66;
	private final TLongObjectHashMap<T> kmers = new TLongObjectHashMap<T>();
	private final int k;
	private DeBruijnGraphExporter<T> exporter = null;
	public DeBruijnGraphBase(int k) {
		this.k = k;
	}
	public int getK() { return k; }
	public T getKmer(long kmer) { return kmers.get(kmer); }
	public int getKmerCount() { return kmers.size(); }
	public Collection<Long> getAllKmers() { return Longs.asList(kmers.keys()); }
	/**
	 * Gets the underlying data structure storing the k-mer graph
	 * @return
	 */
	protected TLongObjectHashMap<T> getBackingStore() { return kmers; }
	public int size() { return kmers.size(); }
	/**
	 * Merges the given nodes together
	 * @param node first node to merge
	 * @param toAdd second node to merge
	 * @return Merged node. This may be a mutable version of one of the arguments 
	 */
	protected abstract T merge(T node, T toAdd);
	/**
	 * Removes the given evidence from a node
	 * @param node node to remove from
	 * @param toAdd second node to merge
	 * @return Merged node. As nodes may be mutable, 
	 */
	protected abstract T remove(T node, T toRemove);
	public T add(T node) {
		T existing = kmers.get(node.kmer());
		if (existing != null) {
			node = merge(existing, node);
			kmers.put(node.kmer(), node);
		} else {
			kmers.put(node.kmer(), node);
			onKmerAdded(node.kmer(), node);
		}
		if (getGraphExporter() != null) {
			getGraphExporter().trackChanges(node.kmer(), node);
		}
		return node;
	}
	/**
	 * removes the given kmer from the graph
	 * @param kmer kmer to remove
	 */
	public void remove(long kmer) {
		kmers.remove(kmer);
		onKmerRemoved(kmer);
	}
	public T remove(long kmer, T node) {
		T existing = kmers.get(kmer);
		if (existing != null) {
			node = remove(existing, node);
			if (node.weight() <= 0) {
				remove(kmer);
			} else {
				kmers.put(kmer, node);
			}
		}
		if (getGraphExporter() != null) {
			getGraphExporter().trackChanges(kmer, node);
		}
		return existing;
	}
	protected void onKmerAdded(long kmer, T node) { }
	protected void onKmerRemoved(long kmer) { }
	/**
	 * Adjusts base qualities to be within valid FASTQ encoding range 
	 * @param bases base qualities to adjust
	 * @return 0-based phred-encodable base qualities
	 */
	public static byte toPicardFastqBaseQuality(int qual) {
		return (byte) (qual > MAX_QUAL_SCORE ? MAX_QUAL_SCORE : qual);
	}
	/**
	 * Gets the all graph kmer following the given kmer
	 * @param state kmer
	 * @param inclusionSet kmers must be in this set. Parameter is ignored if null
	 * @param exclusionSet kmers must not be in this set. Parameter is ignored if null
	 * @return next kmers 
	 */
	public List<Long> nextStates(long state, Set<Long> inclusionSet, Set<Long> exclusionSet) {
		List<Long> result = Lists.newArrayListWithCapacity(4);
		for (Long next : KmerEncodingHelper.nextStates(k, state)) {
			T node = kmers.get(next);
			if (node != null) {
				if ((inclusionSet == null || inclusionSet.contains(next)) && 
						(exclusionSet == null || !exclusionSet.contains(next))) {
					result.add(next);
				}
			}
		}
		return result; 
	}
	/**
	 * Gets the all graph kmer preceeding the given kmer
	 * @param state kmer
	 * @param inclusionSet kmers must be in this set. Parameter is ignored if null
	 * @param exclusionSet kmers must not be in this set. Parameter is ignored if null
	 * @return previous kmers 
	 */
	public List<Long> prevStates(long state, Set<Long> inclusionSet, Set<Long> exclusionSet) {
		List<Long> result = Lists.newArrayListWithCapacity(4);
		for (Long next : KmerEncodingHelper.prevStates(k, state)) {
			T node = kmers.get(next);
			if (node != null) {
				if ((inclusionSet == null || inclusionSet.contains(next)) && 
						(exclusionSet == null || !exclusionSet.contains(next))) {
					result.add(next);
				}
			}
		}
		return result; 
	}
	public Iterable<Long> adjacentStates(long state, Set<Long> inclusionSet, Set<Long> exclusionSet) {
		return Iterables.concat(nextStates(state, inclusionSet, exclusionSet), prevStates(state, inclusionSet, exclusionSet));
	}
	/**
	 * Gets the best kmer following the given kmer
	 * @param state kmer
	 * @param inclusionSet kmers must be in this set. Parameter is ignored if null
	 * @param exclusionSet kmers must not be in this set. Parameter is ignored if null
	 * @return next kmer, null if no valid kmer 
	 */
	public Long greedyNextState(long state, Set<Long> inclusionSet, Set<Long> exclusionSet) {
		int best = Integer.MIN_VALUE;
		Long bestNode = null;
		for (Long next : nextStates(state, inclusionSet, exclusionSet)) {
			int weight = kmers.get(next).weight();
			if (weight > best) {
				bestNode = next;
				best = weight;
			}
		}
		return bestNode; 
	}
	/**
	 * Gets the best kmer preceeding the given kmer
	 * @param state kmer
	 * @param inclusionSet kmers must be in this set. Parameter is ignored if null
	 * @param exclusionSet kmers must not be in this set. Parameter is ignored if null
	 * @return previous kmer, null if no valid kmer 
	 */
	public Long greedyPrevState(long state, Set<Long> inclusionSet, Set<Long> exclusionSet) {
		int best = Integer.MIN_VALUE;
		Long bestNode = null;
		for (Long prev : prevStates(state, inclusionSet, exclusionSet)) {
			int weight = kmers.get(prev).weight();
			if (weight > best) {
				bestNode = prev;
				best = weight;
			}
		}
		return bestNode; 
	}
	/**
	 * Finds the simple path inclusive of the given kmer by traversing to 
	 * the highest weighted adjacent node
	 * @param start seed kmer
	 * @return path
	 */
	public LinkedList<Long> greedyTraverse(Long start) {
		LinkedList<Long> path = new LinkedList<Long>();
		Set<Long> visited = new HashSet<Long>();
		path.add(start);
		visited.add(start);
		for (Long node = greedyPrevState(start, null, visited); node != null; node = greedyPrevState(node, null, visited)) {
			path.addFirst(node);
			visited.add(node);
		}
		for (Long node = greedyNextState(start, null, visited); node != null; node = greedyNextState(node, null, visited)) {
			path.addLast(node);
			visited.add(node);
		}
		return path;
	}
	/**
	 * Base qualities of contig
	 * @param pathAllKmers kmer contig
	 * @return base qualities of a positive strand SAMRecord readout of contig
	 */
	public byte[] getBaseQuals(List<List<T>> pathAllKmers) {
		int[] kmerWeight = new int[pathAllKmers.size()];
		for (int i = 0; i < pathAllKmers.size(); i++) {
			for (T node : pathAllKmers.get(i)) {
				// subtract # reads to adjust for the +1 qual introduced by ReadKmerIterable
				// to ensure positive node weights
				kmerWeight[i] = node.weight() - node.getSupportingEvidenceList().size();
			}
		}
		return kmerWeightsToBaseQuals(k, kmerWeight);
	}
	public static byte[] kmerWeightsToBaseQuals(int k, int[] kmerWeight) {
		byte[] result = new byte[kmerWeight.length + k - 1];
		// take the average kmer qual of the kmers this bases is part of
		int accum = 0;
		int windowSize = 0;
		for (int i = 0; i < result.length; i++) {
			if (i < kmerWeight.length) {
				// add to window
				windowSize++;
				accum += kmerWeight[i];
			}
			if (i >= k) {
				// remove from window
				windowSize--;
				accum -= kmerWeight[i - k];
			}
			result[i] = toPicardFastqBaseQuality(accum / windowSize);
		}
		return result;
	}
	public Set<String> getSupportingEvidence(Iterable<Long> path) {
		Set<String> reads = Sets.newHashSet();
		for (Long kmer : path) {
			reads.addAll(kmers.get(kmer).getSupportingEvidenceList());
		}
		return reads;
	}
	/**
	 * Total kmer support for the given path by category
	 * @param breakendPathAllKmers kmer contig
	 * @return number of read bases per category include in at least one kmer on the given kmer contig
	 */
	public int[] getEvidenceKmerCount(List<List<T>> breakendPathAllKmers) {
		int[] count = new int[1];
		for (T node : Iterables.concat(breakendPathAllKmers)) {
			int[] c = node.getCountByCategory();
			if (c.length > count.length) {
				count = Arrays.copyOf(count, c.length);
			}
			for (int i = 0; i < count.length; i++) {
				count[i] += c[i];
			}
		}
		return count;
	}
	/**
	 * Set of all kmers reachable from the given kmer
	 * @param seed start kmer
	 * @return all reachable kmers
	 */
	public Iterable<Long> reachableFrom(long seed) {
		Set<Long> reachable = Sets.newHashSet();
		Set<Long> frontier = Sets.newHashSet();
		frontier.add(seed);
		while (!frontier.isEmpty()) {
			long kmer = frontier.iterator().next();
			frontier.remove(kmer);
			reachable.add(kmer);
			// Add neighbours of this kmer to the frontier
			for (long adjKmer : adjacentStates(kmer, null, reachable)) {
				// don't need to check for inclusion since frontier is a set
				frontier.add(adjKmer);
			}
		}
		return reachable;
	}
	/**
	 * Ordering of kmers by kmer weight. 
	 */
	public Ordering<Long> ByKmerWeight = new Ordering<Long>() {
		public int compare(Long kmer1, Long kmer2) {
			return Ints.compare(getWeight(kmer1), getWeight(kmer2));
		}
		private int getWeight(long kmer) {
			int weight = 0;
			T node = kmers.get(kmer);
			if (node != null) weight = node.weight();
			return weight;
		}
	};
	@Override
	public String toString() {
		return toString(16);
	}
	public String toString(int maxNodesToPrint) {
		StringBuilder sb = new StringBuilder();
		sb.append(String.format("De Bruijn graph: k=%d, %d kmers\n", k, kmers.size()));
		for (long x : kmers.keys()) {
			sb.append(printKmer(x));
			maxNodesToPrint--;
			if (maxNodesToPrint <= 0) break;
		}
		return sb.toString();
	}
	public String printKmer(long x) {
		StringBuilder sb = new StringBuilder();
		sb.append(String.format("%s(%d): %s",
				KmerEncodingHelper.toString(k, x),
				x,
				kmers.get(x)
				));
		sb.append(" from:{");
		for (Long y : KmerEncodingHelper.prevStates(k, x)) {
			T node = kmers.get(y);
			if (node != null) {
				sb.append(KmerEncodingHelper.toString(k, y));
				sb.append(',');
			}
		}
		sb.append("} to:{");
		for (Long y : KmerEncodingHelper.nextStates(k, x)) {
			T node = kmers.get(y);
			if (node != null) {
				sb.append(KmerEncodingHelper.toString(k, y));
				sb.append(',');
			}
		}
		sb.append("}\n");
		return sb.toString();
	}
	public String debugPrintPaths() {
		Map<Long, Integer> contigLookup = Maps.newHashMap();
		TLongSet remaining = new TLongHashSet(kmers.keySet());
		List<LinkedList<Long>> paths = Lists.newArrayList();
		int contig = 0;
		// enumerate the paths
		while (!remaining.isEmpty()) {
			contig++;
			LinkedList<Long> path = new LinkedList<Long>();
			path.add(remaining.iterator().next());
			remaining.remove(path.getFirst());
			contigLookup.put(path.getFirst(), contig);
			for(List<Long> adj = nextStates(path.getLast(), null, null); adj.size() == 1 && prevStates(adj.get(0), null, null).size() <= 1; adj = nextStates(path.getLast(), null, null)) {
				contigLookup.put(adj.get(0), contig);
				path.addLast(adj.get(0));
				remaining.remove(adj.get(0));
			}
			for(List<Long> adj = prevStates(path.getFirst(), null, null); adj.size() == 1 && nextStates(adj.get(0), null, null).size() <= 1; adj = prevStates(path.getFirst(), null, null)) {
				contigLookup.put(adj.get(0), contig);
				path.addFirst(adj.get(0));
				remaining.remove(adj.get(0));
			}
			paths.add(path);
		}
		return debugPrintPaths(paths, contigLookup);
	}
	private String debugPrintPaths(List<LinkedList<Long>> paths, Map<Long, Integer> contigLookup) {
		StringBuilder sb = new StringBuilder();
		for (LinkedList<Long> path : paths) {
			sb.append(printPath(path, contigLookup));
		}
		return sb.toString();
	}
	private String printPathAttributes(LinkedList<Long> path) {
		return "";
	}
	private String printPath(LinkedList<Long> path, Map<Long, Integer> contigLookup) {
		StringBuilder sb = new StringBuilder();
		sb.append(String.format("[%3d]\t", contigLookup.get(path.getFirst())));
		sb.append(new String(KmerEncodingHelper.baseCalls(path, getK())));
		int weight = 0;
		for (Long y : path) weight += kmers.get(y).weight();
		sb.append(String.format(" %d kmers %d weight", path.size(), weight));
		sb.append(" from:{");
		for (Long y : prevStates(path.getFirst(), null, null)) {
			sb.append(contigLookup.get(y));
			sb.append(',');
		}
		sb.append("} to:{");
		for (Long y : nextStates(path.getLast(), null, null)) {
			sb.append(contigLookup.get(y));
			sb.append(',');
		}
		sb.append("}");
		sb.append(printPathAttributes(path));
		sb.append('\n');
		return sb.toString();
	}
	/**
	 * Gets the graph debugging exporter  
	 * @return
	 */
	public DeBruijnGraphExporter<T> getGraphExporter() {
		return exporter;
	}
	/**
	 * Sets the graph exporter to use for debugging 
	 * @param gexf
	 */
	public void setGraphExporter(DeBruijnGraphExporter<T> exporter) {
		this.exporter = exporter;
	}
	@Override
	public List<T> next(T node) {
		ArrayList<T> result = new ArrayList<T>(2);
		for (long kmer : KmerEncodingHelper.nextStates(k, node.kmer())) {
			T n = kmers.get(kmer);
			if (n != null) {
				result.add(n);
			}
		}
		return result;
	}
	@Override
	public List<T> prev(T node) {
		ArrayList<T> result = new ArrayList<T>(2);
		for (long kmer : KmerEncodingHelper.prevStates(k, node.kmer())) {
			T n = kmers.get(kmer);
			if (n != null) {
				result.add(n);
			}
		}
		return result;
	}
	@Override
	public void removeNode(T node) {
		kmers.remove(node.kmer());
	}
	@Override
	public void removeEdge(T source, T sink) {
		throw new IllegalStateException("Edges are stored implicitly");
	}
	@Override
	public void addNode(T node) {
		add(node);
	}
	@Override
	public void addEdge(T source, T sink) {
		throw new IllegalStateException("Edges are stored implicitly");
	}
	@Override
	public Collection<T> allNodes() {
		return kmers.valueCollection();
	}
	@Override
	public String toString(Iterable<? extends T> path) {
		Iterable<Long> asKmers = Iterables.transform(path, new Function<T, Long>() {
			@Override
			public Long apply(T input) {
				return input.kmer();
			}
		});
		return new String(KmerEncodingHelper.baseCalls(Lists.newArrayList(asKmers), getK()), StandardCharsets.US_ASCII);
	}
	@Override
	public long getKmer(T node) {
		return node.kmer();
	}
	@Override
	public boolean isReference(T node) {
		return node.isReference();
	}
	@Override
	public int getWeight(T node) {
		return node.weight();
	}
}
