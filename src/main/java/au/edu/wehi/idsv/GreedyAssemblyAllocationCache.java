package au.edu.wehi.idsv;

import java.io.IOException;
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
	private final GreedyAllocationCacheLookup<AlignmentScoreNode> bestReadAlignment;
	private final AtomicLong loaded = new AtomicLong(0);
	/**
	 * Creates a new allocation caches
	 * @param threads number of concurrent access threads. This implementation is thread safe for values greater than 1.
	 */
	public GreedyAssemblyAllocationCache(long uniqueReads) {
		bestReadAlignment = createLookup("bestReadAlignment", AlignmentScoreNode.class, uniqueReads);
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
		Hash96bit readkey = new Hash96bit(read);
		Hash96bit alignmentkey = new Hash96bit(alignment);
		putAlignmentScoreNode(bestReadAlignment, readkey, alignmentkey, assemblyScore);
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
		Hash96bit readkey = new Hash96bit(read);
		Hash96bit alignmentkey = new Hash96bit(alignment);
		
		AlignmentScoreNode node = bestReadAlignment.get(readkey);
		return node != null && node.getAlignment().equals(alignmentkey);
	}
	@Override
	public void close() throws IOException {
		bestReadAlignment.close();
	}
}
