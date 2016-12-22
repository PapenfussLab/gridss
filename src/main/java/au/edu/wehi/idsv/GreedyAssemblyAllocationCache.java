package au.edu.wehi.idsv;

import java.util.HashMap;
import java.util.Map.Entry;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.SAMRecord;

/**
 * Tracks the allocation of reads to assemblies to allow
 * allocation of each read to the alignment location resulting in the
 * best assembly. 
 *
 */
public class GreedyAssemblyAllocationCache extends GreedyAllocationCache {
	/**
	 * read -> (event, score, read alignment)
	 * 
	 * Only evidence from the best alignment for each read should be allocated
	 * since we don't want to allocate the mutually exclusive evidence that
	 * results from two separate read alignments for the same read
	 */
	private final HashMap<Hash128bit, Node> bestReadAlignment = new HashMap<>();
	/**
	 * Merges the given caches together 
	 * @param caches caches representing subsets of the assemblies
	 * @return A cache incorporating the full evidence allocation
	 */
	public static GreedyAssemblyAllocationCache merge(Iterable<GreedyAssemblyAllocationCache> caches) {
		GreedyAssemblyAllocationCache merged = new GreedyAssemblyAllocationCache();
		for (GreedyAssemblyAllocationCache cache : caches) {
			for (Entry<Hash128bit, Node> e : cache.bestReadAlignment.entrySet()) {
				Hash128bit key = e.getKey();
				Node value = e.getValue();
				Node existingValue = merged.bestReadAlignment.get(key);
				if (existingValue != null && existingValue.score < value.score) {
					merged.bestReadAlignment.put(key, value);
				}
			}
		}
		return merged;
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
		put(bestReadAlignment, readkey, alignmentkey, null, assemblyScore);
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
		
		Node node = bestReadAlignment.get(readkey);
		return node != null && node.alignment == alignmentkey;
	}
}
