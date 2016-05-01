package au.edu.wehi.idsv.debruijn.positional;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.NonReferenceReadPair;
import au.edu.wehi.idsv.SoftClipEvidence;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.debruijn.PackedKmerList;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.sam.CigarUtil;
import au.edu.wehi.idsv.sam.SAMRecordUtil;

/**
 * Minimal evidence for incorporating soft clip and
 * read pair evidence into a de Bruijn graph
 * 
 * @author Daniel Cameron
 *
 */
public class KmerEvidence extends PackedKmerList {
	// TODO: turn into subclasses
	// - soft clip can store end and BreakendSummary implicitly
	// - read pair can store anchor position & direction instead of breakend & anchor
	// - read pair anchor can point to rp
	private final String id;
	private final BitSet anchor;
	private final BitSet ambiguous;
	private final int start;
	private final int end;
	private final BreakendSummary be;
	private final float score;
	public KmerSupportNode node(int offset) {
		if (ambiguous != null && ambiguous.get(offset)) {
			return null;
		}
		return new KmerSupportNode(this, offset);
	}
	public float evidenceQuality() { return score; }
	public BreakendSummary breakend() { return be; }
	public String evidenceId() { return id; }
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
		if (anchor == null) return false;
		return anchor.get(offset);
	}
	public boolean isAnchored() {
		if (anchor == null) return false;
		return !anchor.isEmpty();
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
	public KmerEvidence(
			String evidenceId,
			int start,
			int end,
			int k, BitSet anchorKmers, byte[] bases, byte[] qual, boolean reverse, boolean complement,
			BreakendSummary be, float evidenceQual) {
		super(k, bases, qual, reverse, complement);
		assert(evidenceId != null);
		assert(qual.length == bases.length);
		this.id = evidenceId;
		this.start = start;
		this.end = end;
		this.ambiguous = ambiguousKmers(k, bases);
		this.anchor = anchorKmers;
		this.be = be;
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
		return new KmerEvidence(pair.getEvidenceID(), startPosition, endPosition, k, null, remote.getReadBases(), remote.getBaseQualities(), reverseComp, reverseComp, pair.getBreakendSummary(), pair.getBreakendQual());
	}
	/**
	 * Creates anchoring evidence for the given read pair
	 * @param k kmer size
	 * @param pair read pair evidence
	 * @return anchoring support
	 */
	public static KmerEvidence createAnchor(int k, NonReferenceReadPair pair, int disallowMismatch, ReferenceLookup reference) {
		return createAnchor(pair.getEvidenceID(), k, pair.getLocalledMappedRead(), pair.getBreakendSummary().direction, disallowMismatch, reference);
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
	public static KmerEvidence createAnchor(String evidenceID, int k, SAMRecord read, BreakendDirection direction, int disallowMismatch, ReferenceLookup reference) {
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
		return new KmerEvidence(evidenceID, firstBasePosition, firstBasePosition, k, anchors, bases, read.getBaseQualities(), false, false, null, 0);
	}
	private static byte getBase(ReferenceLookup reference, int referenceIndex, int contigLength, int position) {
		if (position <= 0 || position > contigLength) return 'N';
		return reference.getBase(referenceIndex, position);
	}
	public static KmerEvidence create(int k, SoftClipEvidence softClipEvidence, boolean trimOtherSoftClip) {
		SAMRecord read = softClipEvidence.getSAMRecord();
		List<CigarElement> elements = read.getCigar().getCigarElements();
		int startClipLength = SAMRecordUtil.getStartSoftClipLength(elements);
		int endClipLength = SAMRecordUtil.getEndSoftClipLength(elements);
		
		if (trimOtherSoftClip && endClipLength > 0 && startClipLength > 0) {
			SAMRecord trimmed;
			try {
				trimmed = (SAMRecord)read.clone();
			} catch (CloneNotSupportedException e) {
				throw new RuntimeException(e);
			}
			List<CigarElement> ce = new ArrayList<CigarElement>(trimmed.getCigar().getCigarElements());
			int trimStart;
			int trimEnd;
			if (softClipEvidence.getBreakendSummary().direction == BreakendDirection.Forward) {
				trimEnd = 0;
				trimStart = startClipLength;
				while (ce.get(0).getOperator() == CigarOperator.SOFT_CLIP || ce.get(0).getOperator() == CigarOperator.HARD_CLIP) {
					ce.remove(0);
				}
				startClipLength = 0;
			} else {
				trimEnd = endClipLength;
				trimStart = 0;
				while (ce.get(ce.size() - 1).getOperator() == CigarOperator.SOFT_CLIP || ce.get(ce.size() - 1).getOperator() == CigarOperator.HARD_CLIP) {
					ce.remove(ce.size() - 1);
				}
				endClipLength = 0;
			}
			trimmed.setCigar(new Cigar(ce));
			trimmed.setBaseQualities(Arrays.copyOfRange(trimmed.getBaseQualities(), trimStart, trimmed.getReadLength() - trimEnd));
			trimmed.setReadBases(Arrays.copyOfRange(trimmed.getReadBases(), trimStart, trimmed.getReadLength() - trimEnd));
			read = trimmed;
		}
		if (k > read.getReadLength()) {
			return null;
		}
		// TODO: handle indels
		BitSet anchor = new BitSet(read.getReadLength());
		int startAnchorIndex = startClipLength;
		int endAnchorIndex = read.getReadLength() - (k - 1) - endClipLength;
		if (endAnchorIndex > startAnchorIndex) {
			anchor.set(startClipLength, read.getReadLength() - (k - 1) - endClipLength, true);
		}
		int startPosition;
		if (softClipEvidence.getBreakendSummary().direction == BreakendDirection.Forward) {
			startPosition = read.getAlignmentEnd() + endClipLength - read.getReadLength() + 1;
		} else {
			startPosition = read.getAlignmentStart() - startClipLength;
		}
		return new KmerEvidence(softClipEvidence.getEvidenceID(), startPosition, startPosition, k, anchor, read.getReadBases(), read.getBaseQualities(), false, false,
				softClipEvidence.getBreakendSummary(), softClipEvidence.getBreakendQual());
	}
	@Override
	public String toString() {
		return evidenceId();
	}
	@Override
	public int hashCode() {
		return id.hashCode();
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
		if (id == null) {
			if (other.id != null)
				return false;
		} else if (!id.equals(other.id))
			return false;
		return true;
	}
}
