package au.edu.wehi.idsv.debruijn.positional;

import htsjdk.samtools.util.Log;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.PriorityQueue;

import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.NonReferenceReadPair;
import au.edu.wehi.idsv.SoftClipEvidence;

import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

/**
 * Transforms a breakend start DirectedEvidence iterator
 * into a start position sorted kmer support iterator
 * @author Daniel Cameron
 *
 */
public class SupportNodeIterator implements PeekingIterator<KmerSupportNode> {
	private static final Log log = Log.getInstance(SupportNodeIterator.class);
	private final PeekingIterator<DirectedEvidence> underlying;
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
	private final PriorityQueue<KmerSupportNode> buffer = new PriorityQueue<KmerSupportNode>(1024, KmerNodeUtil.ByFirstStart);
	private final EvidenceTracker tracker;
	private int inputPosition = Integer.MIN_VALUE;
	private int firstReferenceIndex;
	private int lastPosition = Integer.MIN_VALUE;
	private long consumed = 0;
	private DirectedEvidence lastEvidence = null;
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
	public SupportNodeIterator(int k, Iterator<DirectedEvidence> it, int maxFragmentSize, EvidenceTracker tracker) {
		this.underlying = Iterators.peekingIterator(it);
		this.k = k;
		this.maxSupportStartPositionOffset = maxFragmentSize;
		this.emitOffset = maxSupportStartPositionOffset + 1;
		if (underlying.hasNext()) {
			firstReferenceIndex = underlying.peek().getBreakendSummary().referenceIndex;
		}
		this.tracker = tracker;
	}
	private void process(DirectedEvidence de) {
		if (tracker != null && tracker.isTracked(de.getEvidenceID())) {
			log.error(String.format("Attempting to add %s to assembly when already present. "
					+ "Possible causes are: duplicate read name, alignment with bwa mem without -M flag. "
					+ "Please ensure read names are unique across all inpuut files and there exists only one primary alignment record per read.",
					de.getEvidenceID()));
			return;
		}
		assert(de.getBreakendSummary().referenceIndex == firstReferenceIndex);
		assert(lastEvidence == null || DirectedEvidence.ByStartEnd.compare(lastEvidence, de) <= 0);
		lastEvidence = de;
		
		KmerEvidence e;
		if (de instanceof SoftClipEvidence) {
			e = KmerEvidence.create(k, (SoftClipEvidence)de, true);
		} else if (de instanceof NonReferenceReadPair) {
			e = KmerEvidence.create(k, (NonReferenceReadPair)de);
		} else {
			throw new RuntimeException("Assembler able to process only soft clip and read pair evidence");
		}
		if (e == null) {
			return;
		}
		List<KmerSupportNode> supportNodes = new ArrayList<KmerSupportNode>(e.length());
		boolean hasNonReference = false;
		for (int i = 0; i < e.length(); i++) {
			KmerSupportNode support = e.node(i); 
			if (support != null) {
				// max sure that we are actually able to resort into kmer order
				assert(support.firstStart() >= de.getBreakendSummary().start - maxSupportStartPositionOffset);
				assert(support.weight() > 0);
				supportNodes.add(support);
				hasNonReference |= !support.isReference();
			}
		}
		if (hasNonReference) {
			// only add evidence that proves support for an SV
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
		}
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
			inputPosition = underlying.peek().getBreakendSummary().start;
			advance();
		}
		if (!underlying.hasNext()) {
			inputPosition = Integer.MAX_VALUE;
			advance();
		}
	}
	private void advance() {
		while (underlying.hasNext() && underlying.peek().getBreakendSummary().start <= inputPosition) {
			process(underlying.next());
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