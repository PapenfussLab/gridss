package au.edu.wehi.idsv.debruijn.positional;

import au.edu.wehi.idsv.*;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.debruijn.PackedKmerList;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.sam.CigarUtil;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.util.IntervalUtil;
import au.edu.wehi.idsv.util.MessageThrottler;
import gridss.ComputeSamTags;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.NotImplementedException;
import org.apache.commons.lang3.tuple.Pair;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;

/**
 * Minimal evidence for incorporating soft clip and
 * read pair evidence into a positional de Bruijn graph
 * 
 * @author Daniel Cameron
 *
 */
public class KmerEvidence extends PackedKmerList {
	private static final Log log = Log.getInstance(KmerEvidence.class);
	// TODO: turn into subclasses
	// - soft clip can store end and BreakendSummary implicitly
	// - read pair can store anchor position & direction instead of breakend & anchor
	// - read pair anchor can point to rp
	private final DirectedEvidence evidence;
	private final int refContigLength;
	private final int firstAnchorKmer;
	private final int lastAnchorKmer;
	private BitSet ambiguous;
	private final int start;
	private final int end;
	private final float score;
	private final boolean isReadPairAnchorRead;
	public KmerSupportNode node(int offset) {
		if (ambiguous != null && ambiguous.get(offset)) {
			return null;
		}
		return new KmerSupportNode(this, offset);
	}
	public float evidenceQuality() { return score; }
	public DirectedEvidence evidence() { return evidence; }
	/**
	 * Start position of first kmer
	 * @return
	 */
	public int startPosition() { 
		return start;
	}
	/**
	 * End position of first kmer
	 * @return
	 */
	public int endPosition() { 
		return end;
	}
	public boolean isAnchored(int offset) {
		return offset >= firstAnchorKmer && offset < lastAnchorKmer &&
			start + offset > 0 && end + offset + kmerSize() - 1 <= refContigLength;
	}
	public boolean isAnchored() {
		return firstAnchorKmer < lastAnchorKmer;
	}
	private static BitSet ambiguousKmers(int k, byte[] bases) {
		BitSet ambiguous = null;
		for (int i = 0; i < bases.length; i++) {
			if (KmerEncodingHelper.isAmbiguous(bases[i])) {
				int kmerCount = bases.length - k + 1;
				if (ambiguous == null) {
					ambiguous = new BitSet(kmerCount);
				}
				ambiguous.set(Math.max(0, i - k + 1), Math.min(i, kmerCount - 1) + 1, true);
			}
		}
		return ambiguous;
	}
	private KmerEvidence(
			DirectedEvidence evidence,
			int start,
			int end,
			int k, int firstAnchoredKmer, int lastAnchoredKmer, byte[] bases, byte[] qual, boolean reverse, boolean complement,
			float evidenceQual,
			boolean isReadPairAnchorRead) {
		super(k, bases, qual, reverse, complement);
		assert(evidence != null);
		assert(qual.length == bases.length);
		this.evidence = evidence;
		this.refContigLength = getReferenceContigLength(evidence);
		this.start = start;
		this.end = end;
		this.firstAnchorKmer = firstAnchoredKmer;
		this.lastAnchorKmer = lastAnchoredKmer;
		this.score = evidenceQual;
		this.ambiguous = ambiguousKmers(k, bases);
		if (start != end && evidence.getEvidenceSource().getContext().getConfig().getAssembly().positional.trimSelfIntersectingReads) {
			this.ambiguous = flagSelfIntersectingKmersAsAmbiguous(this.ambiguous);
		}
		this.isReadPairAnchorRead = isReadPairAnchorRead;
	}
	/**
	 * Treats kmers that self-intersect with earlier nodes as ambiguous.
	 * This reduces the explosion nodes that occur in a positional de Bruijn graph
	 * in the presence of low complexity sequence.  
	 */
	private BitSet flagSelfIntersectingKmersAsAmbiguous(BitSet toFlag) {
		// populate lookups
		Long2ObjectOpenHashMap<List<KmerSupportNode>> lookup = new Long2ObjectOpenHashMap<List<KmerSupportNode>>();
		KmerSupportNode[] nodes = new KmerSupportNode[length()];
		for (int i = 0; i < length(); i++) {
			KmerSupportNode n = node(i);
			nodes[i] = n;
			if (n != null) {
				long kmer = n.firstKmer();
				List<KmerSupportNode> kmerList = lookup.get(kmer);
				if (kmerList == null) {
					kmerList = new ArrayList<>(2);
					lookup.put(kmer, kmerList);
				}
				kmerList.add(n);
			}
		}
		// calculcate adjacencies
		List<Pair<Integer, Integer>> unexpectedAdjacencies = new ArrayList<>();
		for (int i = 0; i < length(); i++) {
			KmerSupportNode n = nodes[i];
			if (n != null) {
				long currentkmer = n.firstKmer();
				for (long kmer : KmerEncodingHelper.nextStates(k, currentkmer)) {
					List<KmerSupportNode> kmerList = lookup.get(kmer);
					if (kmerList != null) {
						for (KmerSupportNode adj : kmerList) {
							// we only need to track unexpected adjacencies
							// we already know that it is adjacent to its successor
							if (adj.offset() != n.offset() + 1) {  
								if (IntervalUtil.overlapsClosed(n.firstStart() + 1, n.firstEnd() + 1, adj.firstStart(), adj.firstEnd())) {
									unexpectedAdjacencies.add(Pair.of(n.offset(), adj.offset()));
								}
							}
						}
					}
				}
			}
		}
		// flag the nodes furtherest into the breakpoint as this
		// favours breakpoint truncation over split assemblies
		if (unexpectedAdjacencies.size() > 0 && toFlag == null) {
			toFlag = new BitSet(length());  
		}
		if (evidence.getBreakendSummary().direction == BreakendDirection.Forward) {
			for (Pair<Integer, Integer> pair : unexpectedAdjacencies) {
				toFlag.set(Math.max(pair.getLeft(), pair.getRight()));
			}
		} else {
			for (Pair<Integer, Integer> pair : unexpectedAdjacencies) {
				toFlag.set(Math.min(pair.getLeft(), pair.getRight()));
			}
		}
		return toFlag;
	}
	public static KmerEvidence create(int k, NonReferenceReadPair pair) {
		SAMRecord local = pair.getLocalledMappedRead();
		SAMRecord remote = pair.getNonReferenceRead();
		if (k > remote.getReadLength()) {
			return null;
		}
		boolean reverseComp = !pair.onExpectedStrand();
		int startPosition;
		int endPosition;
		int maxFragSize = pair.getEvidenceSource().getMaxConcordantFragmentSize();
		int minFragSize = pair.getEvidenceSource().getMinConcordantFragmentSize();
		if (pair.getBreakendSummary().direction == BreakendDirection.Forward) {
			//  ----->          anchored
			//  |-------------| max size
			//           <-----
			//  |-----|         minsize
			//   <-----
			int inferredFragStartBasedOnAnchor = local.getAlignmentEnd() - local.getReadLength() + 1 + SAMRecordUtil.getEndSoftClipLength(local);
			int closetPosForFragmentToNotFullyOverlap = local.getAlignmentEnd() - remote.getReadLength() + 2;
			startPosition = Math.max(inferredFragStartBasedOnAnchor + minFragSize - remote.getReadLength(), closetPosForFragmentToNotFullyOverlap);
			endPosition = inferredFragStartBasedOnAnchor + maxFragSize - remote.getReadLength();
		} else {
			//           <----- anchored
			//  |-------------| max size
			//  ----->
			//          |-----| minsize
			//          ----->
			int inferredFragEndBasedOnAnchor = local.getAlignmentStart() + local.getReadLength() - 1 - SAMRecordUtil.getStartSoftClipLength(local);
			int closetPosForFragmentToNotFullyOverlap = local.getAlignmentStart() - 1;
			startPosition = inferredFragEndBasedOnAnchor - maxFragSize + 1;
			endPosition = Math.min(inferredFragEndBasedOnAnchor - minFragSize + 1, closetPosForFragmentToNotFullyOverlap);
		}
		if (remote.getReadBases() == null || remote.getReadBases().length == 0) {
			String msg = String.format("Read %s at %s:%d is missing R2 attribute containing mate information required by GRIDSS. Unable to assemble read",
					local.getReadName(), local.getReferenceName(), local.getAlignmentStart());
			log.error(msg);
			return null;
		}
		if (endPosition < startPosition) {
			if (!MessageThrottler.Current.shouldSupress(log, "kmer filtered")) {
				String msg = String.format("Paired read at %s:%d provides no read pair support. " +
								"The most likely cause is your library fragment size is smaller than your read length. " +
								"Aim for a median library fragment size at least three times your read length when performing library prep.",
						local.getReadName(), local.getReferenceName(), local.getAlignmentStart());
				log.warn(msg);
			}
			return null;
		}
		if (remote.getReadBases() == null || remote.getReadBases().length == 0) {
			String msg = String.format("Read %s at %s:%d is missing R2 attribute containing mate information required by GRIDSS. Unable to assemble read",
					local.getReadName(), local.getReferenceName(), local.getAlignmentStart());
			log.error(msg);
			return null;
		}
		return new KmerEvidence(pair, startPosition, endPosition, k, -1, -1, remote.getReadBases(), remote.getBaseQualities(), reverseComp, reverseComp, pair.getBreakendQual(), false);
	}
	/**
	 * Finds the length of the reference sequence on which this kmer is placed
	 * @param e
	 * @return
	 */
	private static int getReferenceContigLength(DirectedEvidence e) {
		int refIndex = e.getBreakendSummary().referenceIndex;
		EvidenceSource source = e.getEvidenceSource();
		if (source != null) {
			ProcessingContext context = source.getContext();
			if (context != null) {
				return context.getDictionary().getSequence(refIndex).getSequenceLength();
			}
		}
		return Integer.MAX_VALUE;
	}
	/**
	 * Creates anchoring evidence for the given read pair
	 * @param k kmer size
	 * @param pair read pair evidence
	 * @return anchoring support
	 */
	public static KmerEvidence createAnchor(int k, NonReferenceReadPair pair, int disallowMismatch, ReferenceLookup reference) {
		return createAnchor(pair, k, pair.getLocalledMappedRead(), pair.getBreakendSummary().direction, disallowMismatch, reference);
	}
	/**
	 * Creates anchoring evidence for the given read
	 * @param k kmer size
	 * @param read read
	 * @param direction direction to consider matching from. If an indel is present in the read, only bases closes to the
	 * inferred breakend in this direction will be considered anchoring
	 * @param disallowMismatch disallow mismatching bases this number of bases from the end of the read
	 * @param reference reference genome
	 * @return
	 */
	public static KmerEvidence createAnchor(DirectedEvidence evidence, int k, SAMRecord read, BreakendDirection direction, int disallowMismatch, ReferenceLookup reference) {
		if (k > read.getReadLength()) {
			return null;
		}
		BitSet anchors = new BitSet();
		anchors.set(0, read.getReadLength() - k + 1);
		int firstBasePosition;
		read.getCigar().getCigarElements();
		if (direction == BreakendDirection.Forward) {
			firstBasePosition = read.getUnclippedEnd() - read.getReadLength() + 1;
		} else {
			firstBasePosition = read.getUnclippedStart();
		}
		List<CigarElement> cigar = CigarUtil.asUngapped(read.getCigar(), direction == BreakendDirection.Backward);
		byte[] bases = Arrays.copyOf(read.getReadBases(), read.getReadLength());
		// consider clipped bases as ambiguous
		int readOffset = 0;
		for (CigarElement ci : cigar) {
			if (ci.getOperator() == CigarOperator.SOFT_CLIP) {
				for (int i = 0; i < ci.getLength(); i++) {
					bases[readOffset + i] = 'N';
				}
			}
			if (ci.getOperator().consumesReadBases()) {
				readOffset += ci.getLength();
			}
		}
		if (reference != null) {
			int contigLength = reference.getSequenceDictionary().getSequence(read.getReferenceIndex()).getSequenceLength();
			// consider mismatches at end of read ambiguous
			for (int i = 0; i < disallowMismatch && i < bases.length; i++) {
				if (bases[i] != getBase(reference, read.getReferenceIndex(), contigLength, firstBasePosition + i)) {
					bases[i] = 'N';
				}
				int endOffset = bases.length - 1 - i;
				if (bases[endOffset] != getBase(reference, read.getReferenceIndex(), contigLength, firstBasePosition + endOffset)) {
					bases[endOffset] = 'N';
				}
			}
		}
		return new KmerEvidence(evidence, firstBasePosition, firstBasePosition, k, 0, bases.length, bases, read.getBaseQualities(), false, false, 0, true);
	}
	private static byte getBase(ReferenceLookup reference, int referenceIndex, int contigLength, int position) {
		if (position <= 0 || position > contigLength) return 'N';
		return reference.getBase(referenceIndex, position);
	}
	public static KmerEvidence create(int k, SingleReadEvidence sre) {
		if (!sre.isBreakendExact()) {
			throw new NotImplementedException("reassembly of XNX placeholder contigs");
		}
		if (sre.getSAMRecord().getTransientAttribute("HC") != null) {
			if (!MessageThrottler.Current.shouldSupress(log, "hard clipped bases")) {
				log.warn(String.format("Read %s is hard clipped. Please run %s to soften hard clips.",
						sre.getSAMRecord().getReadName(), ComputeSamTags.class.getName()));
			}
		}
		byte[] aseq = sre.getAnchorSequence();
		byte[] beseq = sre.getBreakendSequence();
		byte[] aqual = sre.getAnchorQuality();
		byte[] bequal = sre.getBreakendQuality();
		
		byte[] seq;
		byte[] qual;
		BreakendSummary bs = sre.getBreakendSummary();
		int positionOffset;
		int firstAnchoredBase;
		int anchoredBases = aseq.length;
		if (bs.direction == BreakendDirection.Forward) {
			seq = ArrayUtils.addAll(aseq, beseq);
			qual =  ArrayUtils.addAll(aqual, bequal);
			positionOffset = -(aseq.length - 1);
			firstAnchoredBase = 0;
		} else {
			seq = ArrayUtils.addAll(beseq, aseq);
			qual =  ArrayUtils.addAll(bequal, aqual);
			firstAnchoredBase = beseq.length;
			positionOffset = -beseq.length;
		}
		if (k > seq.length) {
			return null;
		}
		return new KmerEvidence(sre, bs.start + positionOffset, bs.end + positionOffset, k, firstAnchoredBase, firstAnchoredBase + anchoredBases - (k - 1), seq, qual, false, false, sre.getBreakendQual(), false);
		
	}
	@Override
	public String toString() {
		return evidence.getEvidenceID();
	}
	@Override
	public int hashCode() {
		return evidence.getEvidenceID().hashCode() + start + (isReadPairAnchorRead ? 1 : 0);
	}
	public boolean equals(KmerEvidence other) {
		return isReadPairAnchorRead == other.isReadPairAnchorRead &&
				start == other.start &&
				evidence.getEvidenceID().equals(other.evidence.getEvidenceID());
	}
	@Override
	public boolean equals(Object obj) {
		return this == obj || (obj instanceof KmerEvidence && equals((KmerEvidence)obj));
	}
}
