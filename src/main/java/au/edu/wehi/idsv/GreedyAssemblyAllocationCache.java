package au.edu.wehi.idsv;

import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicLong;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;

/**
 * Tracks the allocation of reads to assemblies to allow
 * allocation of each read to the alignment location resulting in the
 * best assembly. 
 *
 */
public class GreedyAssemblyAllocationCache extends GreedyAllocationCache {
	private static final Log log = Log.getInstance(GreedyAssemblyAllocationCache.class);
	/**
	 * read -> (event, score, read alignment)
	 * 
	 * Only evidence from the best alignment for each read should be allocated
	 * since we don't want to allocate the mutually exclusive evidence that
	 * results from two separate read alignments for the same read
	 */
	private final Map<Hash128bit, AlignmentScoreNode> bestReadAlignment;
	private final AtomicLong loaded = new AtomicLong(0);
	public GreedyAssemblyAllocationCache() {
		this(1);
	}
	/**
	 * Creates a new allocation caches
	 * @param threads number of concurrent access threads. This implementation is thread safe for values greater than 1.
	 */
	public GreedyAssemblyAllocationCache(int threads) {
		if (threads <= 1) {
			bestReadAlignment = new HashMap<>(INITIAL_HASH_MAP_SIZE);
		} else {
			bestReadAlignment = new ConcurrentHashMap<Hash128bit, AlignmentScoreNode>(INITIAL_HASH_MAP_SIZE, 0.75f, threads);
		}
	}
	/**
	 * Merges the given caches together 
	 * @param caches caches representing subsets of the assemblies
	 * @return A cache incorporating the full evidence allocation
	 */
	public static GreedyAssemblyAllocationCache merge(Iterable<GreedyAssemblyAllocationCache> caches) {
		GreedyAssemblyAllocationCache merged = new GreedyAssemblyAllocationCache();
		for (GreedyAssemblyAllocationCache cache : caches) {
			merged.addAll(cache);
		}
		return merged;
	}
	public void addAll(GreedyAssemblyAllocationCache cache) {
		for (Entry<Hash128bit, AlignmentScoreNode> e : cache.bestReadAlignment.entrySet()) {
			Hash128bit key = e.getKey();
			AlignmentScoreNode value = e.getValue();
			AlignmentScoreNode existingValue = bestReadAlignment.get(key);
			if (existingValue == null || existingValue.getScore() < value.getScore()) {
				bestReadAlignment.put(key, value);
			}
		}
	}
	protected void addBreakendAssemblyAllocation(float assemblyScore, DirectedEvidence evidence) {
		SAMRecord anchor;
		if (evidence instanceof NonReferenceReadPair) {
			anchor = ((NonReferenceReadPair)evidence).getLocalledMappedRead();
		} else {
			anchor = ((SingleReadEvidence)evidence).getSAMRecord();
		}
		String read = anchor.getReadName() + "#" + Integer.valueOf(SAMRecordUtil.getSegmentIndex(anchor));
		String alignment = getReadAlignment(anchor);
		Hash128bit readkey = new Hash128bit(read);
		Hash128bit alignmentkey = new Hash128bit(alignment);
		put(bestReadAlignment, readkey, alignmentkey, assemblyScore);
		long count = loaded.incrementAndGet();
		if (count % 1000000 == 0) {
			log.info(String.format("Loaded %,d records. Current java heap memory usage is %,d MiB", count, (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) >> 20));
		}
	}
	public boolean isBestBreakendAssemblyAllocation(DirectedEvidence evidence) {
		SAMRecord anchor;
		if (evidence instanceof NonReferenceReadPair) {
			anchor = ((NonReferenceReadPair)evidence).getLocalledMappedRead();
		} else {
			anchor = ((SingleReadEvidence)evidence).getSAMRecord();
		}
		String read = anchor.getReadName() + "#" + Integer.valueOf(SAMRecordUtil.getSegmentIndex(anchor));
		String alignment = getReadAlignment(anchor);
		Hash128bit readkey = new Hash128bit(read);
		Hash128bit alignmentkey = new Hash128bit(alignment);
		
		AlignmentScoreNode node = bestReadAlignment.get(readkey);
		return node != null && node.getAlignment().equals(alignmentkey);
	}
}
