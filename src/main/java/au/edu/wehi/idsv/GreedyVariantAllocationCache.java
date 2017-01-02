package au.edu.wehi.idsv;

import java.util.HashMap;
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
	private final HashMap<Hash128bit, EventAlignmentScoreNode> bestReadPairAlignment;
	/**
	 * read -> (event, score, read alignment)
	 * 
	 * Only evidence from the best alignment for each read should be allocated
	 * since we don't want to allocate the mutually exclusive evidence that
	 * results from two separate read alignments for the same read
	 */
	private final HashMap<Hash128bit, EventAlignmentScoreNode> bestReadAlignment;
	/**
	 * evidenceID -> best event lookup
	 * Each evidence can support only a single variant.
	 */
	private final HashMap<Hash128bit, EventScoreNode> bestEventForEvidence;
	private final AtomicLong loaded = new AtomicLong(0);
	public GreedyVariantAllocationCache(boolean ensureUniqueReadPairAlignment, boolean ensureUniqueReadAlignment, boolean ensureUniqueEvidenceAllocation) {
		this.bestReadPairAlignment = ensureUniqueReadPairAlignment ? new HashMap<>(INITIAL_HASH_MAP_SIZE) : null;
		this.bestReadAlignment = ensureUniqueReadAlignment ? new HashMap<>(INITIAL_HASH_MAP_SIZE) : null;
		this.bestEventForEvidence = ensureUniqueEvidenceAllocation ? new HashMap<>(1000000) : null;
	}
	private static Hash128bit getEvent(VariantContextDirectedBreakpoint variant) {
		return new Hash128bit(variant.getAttributeAsString(VcfSvConstants.BREAKEND_EVENT_ID_KEY, null));
	}
	public void addBreakpoint(String event, float score, DirectedEvidence evidence) {
		addBreakpoint(new Hash128bit(event), score, evidence);
		long count = loaded.incrementAndGet();
		if (count % 1000000 == 0) {
			log.info(String.format("Loaded %,d records. Current java heap memory usage is %,d MiB", count, (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) >> 20));
		}
	}
	public void addBreakpoint(VariantContextDirectedBreakpoint variant, List<DirectedEvidence> evidence) {
		Hash128bit event = getEvent(variant);
		for (DirectedEvidence e : evidence) {
			addBreakpoint(event, variant.getBreakpointQual(), e);
		}
	}
	public boolean isBestBreakpoint(VariantContextDirectedBreakpoint variant, DirectedEvidence evidence) {
		return isBestBreakpoint(getEvent(variant), evidence);
	}
	public boolean isBestBreakpoint(String event, DirectedEvidence evidence) {
		return isBestBreakpoint(new Hash128bit(event), evidence);
	}
	protected void addBreakpoint(Hash128bit event, float score, DirectedEvidence evidence) {
		put(bestEventForEvidence, new Hash128bit(evidence.getEvidenceID()), event, score);
		if (evidence instanceof NonReferenceReadPair) {
			NonReferenceReadPair dp = (NonReferenceReadPair)evidence;
			Hash128bit readpairid = new Hash128bit(dp.getLocalledMappedRead().getReadName());
			Hash128bit alignment = getReadPairAlignment(dp.getLocalledMappedRead());
			put(bestReadPairAlignment, readpairid, alignment, event, score);
		} else {
			assert(evidence instanceof SingleReadEvidence);
			SingleReadEvidence sre = (SingleReadEvidence)evidence;
			SAMRecord r = sre.getSAMRecord();
			Hash128bit readid = new Hash128bit(r.getReadName() + "#" + Integer.valueOf(SAMRecordUtil.getSegmentIndex(r)));
			Hash128bit alignment = new Hash128bit(getReadAlignment(r));
			put(bestReadAlignment, readid, alignment, event, score);
		}
	}
	public boolean isBestBreakpoint(Hash128bit event, DirectedEvidence evidence) {
		if (bestEventForEvidence != null && !event.equals(bestEventForEvidence.get(new Hash128bit(evidence.getEvidenceID())).getEvent())) {
			// This is not the best breakpoint supported by this evidence
			return false;
		}
		if (evidence instanceof NonReferenceReadPair) {
			NonReferenceReadPair dp = (NonReferenceReadPair)evidence;
			Hash128bit readpairid = new Hash128bit(dp.getLocalledMappedRead().getReadName());
			Hash128bit alignment = getReadPairAlignment(dp.getLocalledMappedRead());
			return isBestAlignment(bestReadPairAlignment, readpairid, alignment);
		} else {
			assert(evidence instanceof SingleReadEvidence);
			SingleReadEvidence sre = (SingleReadEvidence)evidence;
			SAMRecord r = sre.getSAMRecord();
			Hash128bit readid = new Hash128bit(r.getReadName() + "#" + Integer.valueOf(SAMRecordUtil.getSegmentIndex(r)));
			Hash128bit alignment = new Hash128bit(getReadAlignment(r));
			return isBestAlignment(bestReadAlignment, readid, alignment);
		}
	}
}
