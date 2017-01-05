package au.edu.wehi.idsv;

import java.io.IOException;
import java.util.List;
import java.util.concurrent.atomic.AtomicLong;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.vcf.VcfSvConstants;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;

public class GreedyVariantAllocationCache extends GreedyAllocationCache {
	private static final Log log = Log.getInstance(GreedyVariantAllocationCache.class);
	/**
	 * read (pair) -> (event, score, read pair alignment)
	 * Only best placement of the read pairs should be allocated.
	 */
	private final GreedyAllocationCacheLookup<EventAlignmentScoreNode> bestReadPairAlignment;
	/**
	 * read -> (event, score, read alignment)
	 * 
	 * Only evidence from the best alignment for each read should be allocated
	 * since we don't want to allocate the mutually exclusive evidence that
	 * results from two separate read alignments for the same read
	 */
	private final GreedyAllocationCacheLookup<EventAlignmentScoreNode> bestReadAlignment;
	/**
	 * evidenceID -> best event lookup
	 * Each evidence can support only a single variant.
	 */
	private final GreedyAllocationCacheLookup<EventScoreNode> bestEventForEvidence;
	private final AtomicLong loaded = new AtomicLong(0);
	public GreedyVariantAllocationCache(boolean ensureUniqueReadPairAlignment, boolean ensureUniqueReadAlignment, boolean ensureUniqueEvidenceAllocation, long uniqueReads) {
		this.bestReadPairAlignment = ensureUniqueReadPairAlignment ? createLookup("bestReadPairAlignment", EventAlignmentScoreNode.class, uniqueReads) : null;
		this.bestReadAlignment = ensureUniqueReadAlignment ? createLookup("bestReadAlignment", EventAlignmentScoreNode.class, uniqueReads) : null;
		this.bestEventForEvidence = ensureUniqueEvidenceAllocation ? createLookup("bestEventForEvidence", EventScoreNode.class, uniqueReads) : null;
	}
	private static Hash96bit getEvent(VariantContextDirectedBreakpoint variant) {
		return new Hash96bit(variant.getAttributeAsString(VcfSvConstants.BREAKEND_EVENT_ID_KEY, null));
	}
	public void addBreakpoint(String event, float score, DirectedEvidence evidence) {
		addBreakpoint(new Hash96bit(event), score, evidence);
		
	}
	public void addBreakpoint(VariantContextDirectedBreakpoint variant, List<DirectedEvidence> evidence) {
		Hash96bit event = getEvent(variant);
		for (DirectedEvidence e : evidence) {
			addBreakpoint(event, variant.getBreakpointQual(), e);
		}
	}
	public boolean isBestBreakpoint(VariantContextDirectedBreakpoint variant, DirectedEvidence evidence) {
		return isBestBreakpoint(getEvent(variant), evidence);
	}
	public boolean isBestBreakpoint(String event, DirectedEvidence evidence) {
		return isBestBreakpoint(new Hash96bit(event), evidence);
	}
	protected void addBreakpoint(Hash96bit event, float score, DirectedEvidence evidence) {
		putEventScoreNode(bestEventForEvidence, new Hash96bit(evidence.getEvidenceID()), event, score);
		if (evidence instanceof NonReferenceReadPair) {
			NonReferenceReadPair dp = (NonReferenceReadPair)evidence;
			Hash96bit readpairid = new Hash96bit(dp.getLocalledMappedRead().getReadName());
			Hash96bit alignment = getReadPairAlignment(dp.getLocalledMappedRead());
			putEventAlignmentScoreNode(bestReadPairAlignment, readpairid, alignment, event, score);
		} else {
			assert(evidence instanceof SingleReadEvidence);
			SingleReadEvidence sre = (SingleReadEvidence)evidence;
			SAMRecord r = sre.getSAMRecord();
			Hash96bit readid = new Hash96bit(r.getReadName() + "#" + Integer.valueOf(SAMRecordUtil.getSegmentIndex(r)));
			Hash96bit alignment = new Hash96bit(getReadAlignment(r));
			putEventAlignmentScoreNode(bestReadAlignment, readid, alignment, event, score);
		}
		long count = loaded.incrementAndGet();
		if (count % 1000000 == 0) {
			log.info(String.format("Loaded %,d records. Current java heap memory usage is %,d MiB",
					count,
					(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) >> 20));
		}
	}
	public boolean isBestBreakpoint(Hash96bit event, DirectedEvidence evidence) {
		if (bestEventForEvidence != null) {
			EventScoreNode lookup = bestEventForEvidence.get(new Hash96bit(evidence.getEvidenceID()));
			if (lookup == null || !event.equals(lookup.getEvent())) {
				// This is not the best breakpoint supported by this evidence
				return false;
			}
		}
		if (evidence instanceof NonReferenceReadPair) {
			NonReferenceReadPair dp = (NonReferenceReadPair)evidence;
			Hash96bit readpairid = new Hash96bit(dp.getLocalledMappedRead().getReadName());
			Hash96bit alignment = getReadPairAlignment(dp.getLocalledMappedRead());
			return isBestAlignment(bestReadPairAlignment, readpairid, alignment);
		} else {
			assert(evidence instanceof SingleReadEvidence);
			SingleReadEvidence sre = (SingleReadEvidence)evidence;
			SAMRecord r = sre.getSAMRecord();
			Hash96bit readid = new Hash96bit(r.getReadName() + "#" + Integer.valueOf(SAMRecordUtil.getSegmentIndex(r)));
			Hash96bit alignment = new Hash96bit(getReadAlignment(r));
			return isBestAlignment(bestReadAlignment, readid, alignment);
		}
	}
	@Override
	public void close() throws IOException {
		bestEventForEvidence.close();
		bestReadAlignment.close();
		bestReadPairAlignment.close();
	}
}
