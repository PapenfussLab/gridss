package au.edu.wehi.idsv.debruijn.positional;

import au.edu.wehi.idsv.Defaults;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.NonReferenceReadPair;
import au.edu.wehi.idsv.SingleReadEvidence;
import au.edu.wehi.idsv.debruijn.positional.optimiseddatastructures.KmerNodeByFirstStartPriorityQueue;
import au.edu.wehi.idsv.util.MessageThrottler;
import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;

import java.util.*;

import static au.edu.wehi.idsv.Defaults.SANITY_CHECK_EVIDENCE_TRACKER;

/**
 * Transforms a breakend start DirectedEvidence iterator
 * into a start position sorted kmer support iterator
 * @author Daniel Cameron
 *
 */
public class SupportNodeIterator implements PeekingIterator<KmerSupportNode> {
	private static final Log log = Log.getInstance(SupportNodeIterator.class);
	private final PeekingIterator<DirectedEvidence> underlying;
	private final boolean includePairAnchors;
	private final int k;
	/**
	 * Position to emit kmers as no more can be added
	 *                  | input position
	 *                  v
	 *             ?????------------<====  RP
	 *              MMMMMSSSSS SC
	 *              SSSSSMMMMM SC
	 *              ====>------------?????  RP
	 */
	private final int emitOffset;
	private final int maxSupportStartPositionOffset;
	private final Queue<KmerSupportNode> buffer = Defaults.USE_OPTIMISED_ASSEMBLY_DATA_STRUCTURES ? new KmerNodeByFirstStartPriorityQueue<>(16) : new PriorityQueue<>(1024, KmerNodeUtil.ByFirstStart);
	private final EvidenceTracker tracker;
	private final int disallowMismatch;
	private int inputPosition = Integer.MIN_VALUE;
	private int firstReferenceIndex;
	private int lastPosition = Integer.MIN_VALUE;
	private long consumed = 0;
	/**
	 * Iterator that converts evidence to kmer nodes 
	 * @param k kmer
	 * @param it underlying evidence iterator sorted by evidence start position
	 * @param maxSupportStartPositionOffset maximum number of bases offset
	 * For typical Illumina PE data, this is equal to the maximum concordant fragment size
	 * This occurs when:
	 *  - short DP causes 1bp breakend interval
	 *  - 1bp of read mapped
	 * In the above scenario, the mate will be placed for assembly purposes as far away
	 * as max frag size from the mapping position (which in this scenario is also the breakend
	 * position).
	 */
	public SupportNodeIterator(int k, Iterator<DirectedEvidence> it, int maxSupportStartPositionOffset, EvidenceTracker tracker, boolean includePairAnchors, int disallowMismatch) {
		this.underlying = Iterators.peekingIterator(it);
		this.k = k;
		this.includePairAnchors = includePairAnchors;
		this.disallowMismatch = disallowMismatch;
		this.maxSupportStartPositionOffset = maxSupportStartPositionOffset;
		this.emitOffset = maxSupportStartPositionOffset + 1 + maxSupportStartPositionOffset; // Need to adjust for evidence in SAMRecord sort order
		if (underlying.hasNext()) {
			firstReferenceIndex = underlying.peek().getBreakendSummary().referenceIndex;
		}
		this.tracker = tracker;
	}
	private void process(DirectedEvidence de) {
		if (tracker != null && tracker.isTracked(de.getEvidenceID())) {
			if (!MessageThrottler.Current.shouldSupress(log, "assembly duplicated reads")) {
				log.warn(String.format("Attempting to add %s (from %s) to assembly when already present. "
						+ "Possible causes are: duplicate read name, alignment with multi-mapping aligner which writes read alignments as distinct pairs. ",
						de.getEvidenceID(), de.getUnderlyingSAMRecord().getReadName()));
			}
			return;
		}
		assert(de.getBreakendSummary().referenceIndex == firstReferenceIndex);
		/*
		if (lastEvidence != null && DirectedEvidence.ByStartEnd.compare(lastEvidence, de) > 0) {
			String msg = String.format("SupportNodeIterator requires evidence to be sorted by starting position. Encountered %s at %s before %s at %s.",
					lastEvidence.getEvidenceID(), lastEvidence.getBreakendSummary(),
					de.getEvidenceID(), de.getBreakendSummary());
			log.error(msg);
			throw new RuntimeException(msg);
		}
		lastEvidence = de;
		*/
		KmerEvidence e;
		KmerEvidence e2 = null;
		if (de instanceof SingleReadEvidence) {
			e = KmerEvidence.create(k, (SingleReadEvidence)de);
		} else if (de instanceof NonReferenceReadPair) {
			NonReferenceReadPair nrrp = (NonReferenceReadPair)de;
			e = KmerEvidence.create(k, nrrp);
			if (includePairAnchors) {
				e2 = KmerEvidence.createAnchor(k, nrrp, disallowMismatch, nrrp.getEvidenceSource().getContext().getReference());
			}
		} else {
			throw new RuntimeException("Assembler able to process only soft clip and read pair evidence");
		}
		if (e == null) {
			return;
		}
		List<KmerSupportNode> supportNodes = new ArrayList<KmerSupportNode>(e.length() + (e2 == null ? 0 : e2.length()));
		boolean hasNonReference = addSupport(supportNodes, de, e);
		addSupport(supportNodes, de, e2);
		if (hasNonReference) {
			// only add evidence that provides support for an SV
			// If we have no non-reference kmers then we might
			// never call a contig containing this evidence thus
			// never remove it from the graph
			// SC or RPs with no non-reference kmers can occur when
			// an ambiguous base case exist in the soft clip/mate  
			buffer.addAll(supportNodes);
			if (tracker != null) {
				for (KmerSupportNode sn : supportNodes) {
					tracker.track(sn);
				}
			}
		} else {
			log.debug("Ref anchor");
		}
		if (SANITY_CHECK_EVIDENCE_TRACKER) {
			tracker.sanityCheck();
		}
	}
	private boolean addSupport(List<KmerSupportNode> supportNodes, DirectedEvidence de, KmerEvidence e) {
		boolean hasNonReference = false;
		if (e != null) {
			for (int i = 0; i < e.length(); i++) {
				KmerSupportNode support = e.node(i); 
				if (support != null) {
					// make sure that we are actually able to resort into kmer order
					if (support.firstStart() < de.getBreakendSummary().start - maxSupportStartPositionOffset) {
						SAMRecord read = null;
						if (de instanceof SingleReadEvidence) {
							read = ((SingleReadEvidence)de).getSAMRecord(); 
						} else if (de instanceof NonReferenceReadPair) {
							read = ((NonReferenceReadPair)de).getLocalledMappedRead();
						}
						String readString = "";
						if (read != null) {
							readString = read.getReadName();
							if (!read.getReadUnmappedFlag()) {
								readString += String.format(" (%s:%d %s)", read.getReferenceName(), read.getStart(), read.getCigarString());
							}
						}
						String msg = String.format("Error: kmer in evidence %s of read %s out of bounds."
								+ " Kmer support starts at %d which is more than %d before the breakpoint start position at %s",
								de.getEvidenceID(),
								readString,
								support.firstStart(), maxSupportStartPositionOffset, de.getBreakendSummary());
						log.error(msg);
						// Try to continue
						//throw new RuntimeException(msg);
					} else if (support.weight() <= 0) {
						String msg = String.format("Invalid support weight of %d for evidence %s", support.weight(), de.getEvidenceID());
						log.error(msg);
						throw new RuntimeException(msg);
					} else {
						assert(support.weight() > 0);
						supportNodes.add(support);
						hasNonReference |= !support.isReference();
					}
				}
			}
		}
		return hasNonReference;
	}
	@Override
	public boolean hasNext() {
		ensureBuffer();
		return !buffer.isEmpty();
	}
	@Override
	public KmerSupportNode next() {
		ensureBuffer();
		KmerSupportNode node = buffer.poll();
		assert(node.lastStart() >= lastPosition);
		lastPosition = node.lastStart();		
		return node;
	}
	@Override
	public KmerSupportNode peek() {
		ensureBuffer();
		return buffer.peek();
	}
	private void ensureBuffer() {
		while (underlying.hasNext() && (buffer.isEmpty() || buffer.peek().lastStart() > inputPosition - emitOffset)) {
			inputPosition = underlying.peek().getUnderlyingSAMRecord().getAlignmentStart();
			advance();
		}
		if (!underlying.hasNext()) {
			inputPosition = Integer.MAX_VALUE;
			advance();
		}
	}
	private void advance() {
		while (underlying.hasNext() && underlying.peek().getUnderlyingSAMRecord().getAlignmentStart() <= inputPosition) {
			DirectedEvidence nextRecord = underlying.next();
			process(nextRecord);
			consumed++;
		}
	}

	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}
	public int tracking_processedSize() {
		return buffer.size();
	}
	public int tracking_inputPosition() {
		return inputPosition;
	}
	public long tracking_underlyingConsumed() {
		return consumed;
	}
}