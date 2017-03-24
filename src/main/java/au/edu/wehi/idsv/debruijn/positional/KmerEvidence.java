package au.edu.wehi.idsv.debruijn.positional;

import java.util.Arrays;
import java.util.BitSet;
import java.util.List;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.NotImplementedException;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.EvidenceSource;
import au.edu.wehi.idsv.NonReferenceReadPair;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SingleReadEvidence;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.debruijn.PackedKmerList;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.sam.CigarUtil;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

/**
 * Minimal evidence for incorporating soft clip and
 * read pair evidence into a positional de Bruijn graph
 * 
 * @author Daniel Cameron
 *
 */
public class KmerEvidence extends PackedKmerList {
	// TODO: turn into subclasses
	// - soft clip can store end and BreakendSummary implicitly
	// - read pair can store anchor position & direction instead of breakend & anchor
	// - read pair anchor can point to rp
	private final DirectedEvidence evidence;
	private final int refContigLength;
	private final int firstAnchorKmer;
	private final int lastAnchorKmer;
	private final BitSet ambiguous;
	private final int start;
	private final int end;
	private final float score;
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
			float evidenceQual) {
		super(k, bases, qual, reverse, complement);
		assert(evidence != null);
		assert(qual.length == bases.length);
		this.evidence = evidence;
		this.refContigLength = getReferenceContigLength(evidence);
		this.start = start;
		this.end = end;
		this.ambiguous = ambiguousKmers(k, bases);
		this.firstAnchorKmer = firstAnchoredKmer;
		this.lastAnchorKmer = lastAnchoredKmer;
		this.score = evidenceQual;
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
			startPosition = local.getUnclippedStart() + minFragSize - remote.getReadLength();
			endPosition = local.getUnclippedStart() + maxFragSize - remote.getReadLength();
		} else {
			//           <----- anchored
			//  |-------------| max size
			//  ----->
			//          |-----| minsize
			//          -----> 
			startPosition = local.getUnclippedEnd() - maxFragSize + 1;
			endPosition = local.getUnclippedEnd() - minFragSize + 1;
		}
		return new KmerEvidence(pair, startPosition, endPosition, k, -1, -1, remote.getReadBases(), remote.getBaseQualities(), reverseComp, reverseComp, pair.getBreakendQual());
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
				return context.getReference().getSequenceDictionary().getSequence(refIndex).getSequenceLength();
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
		return new KmerEvidence(evidence, firstBasePosition, firstBasePosition, k, 0, bases.length, bases, read.getBaseQualities(), false, false, 0);
	}
	private static byte getBase(ReferenceLookup reference, int referenceIndex, int contigLength, int position) {
		if (position <= 0 || position > contigLength) return 'N';
		return reference.getBase(referenceIndex, position);
	}
	public static KmerEvidence create(int k, SingleReadEvidence sre) {
		if (!sre.isBreakendExact()) {
			throw new NotImplementedException("reassembly of XNX placeholder contigs");
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
		return new KmerEvidence(sre, bs.start + positionOffset, bs.end + positionOffset, k, firstAnchoredBase, firstAnchoredBase + anchoredBases - (k - 1), seq, qual, false, false, sre.getBreakendQual());
		
	}
	@Override
	public String toString() {
		return evidence.getEvidenceID();
	}
	@Override
	public int hashCode() {
		return evidence.getEvidenceID().hashCode();
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		KmerEvidence other = (KmerEvidence) obj;
		String id = evidence.getEvidenceID();
		String oid = other.evidence.getEvidenceID();
		if (id == null) {
			if (oid != null)
				return false;
		} else if (!id.equals(oid))
			return false;
		return true;
	}
}
