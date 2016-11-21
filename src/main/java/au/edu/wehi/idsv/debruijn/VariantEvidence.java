package au.edu.wehi.idsv.debruijn;

import java.util.BitSet;
import java.util.List;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.LinearGenomicCoordinate;
import au.edu.wehi.idsv.NonReferenceReadPair;
import au.edu.wehi.idsv.SingleReadEvidence;
import au.edu.wehi.idsv.SoftClipEvidence;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

/**
 * Helper wrapper exposing information on both the raw evidence read and
 * kmer-specific information for that read.
 *   
 * @author Daniel Cameron
 *
 */
public class VariantEvidence {
	private final PackedKmerList kmers;
	private final BitSet ambiguous;
	/**
	 * Expected linear genomic position of the start of the read
	 */
	private final long expectedAnchorPos;
	/**
	 * Offset of first reference kmer (inclusive)
	 */
	private final int referenceKmerStartOffset;
	/**
	 * Offset of last reference kmer (exclusive)
	 */
	private final int referenceKmerEndOffset;
	private final int startSkipKmerCount;
	private final int endSkipKmerCount;
	private final String evidenceID;
	private boolean isExact;
	private BreakendSummary be;
	private int category;
	public VariantEvidence(int k, NonReferenceReadPair pair, LinearGenomicCoordinate lgc) {
		this.evidenceID = pair.getEvidenceID();
		boolean shouldReverseComplement = !pair.onExpectedStrand();		
		byte[] bases = pair.getNonReferenceRead().getReadBases();
		byte[] quals = pair.getNonReferenceRead().getBaseQualities();
		int chrPos;
		if (pair.getBreakendSummary().direction == BreakendDirection.Forward) {
			// Assumes FR orientation
			chrPos = pair.getLocalledMappedRead().getUnclippedStart() + pair.getEvidenceSource().getExpectedFragmentSize() - bases.length;
		} else {
			chrPos = pair.getLocalledMappedRead().getUnclippedEnd() - pair.getEvidenceSource().getExpectedFragmentSize() + 1;
		}
		this.expectedAnchorPos = lgc.getLinearCoordinate(pair.getBreakendSummary().referenceIndex, chrPos);
		this.referenceKmerStartOffset = -1;
		this.referenceKmerEndOffset = -1;
		this.startSkipKmerCount = 0;
		this.endSkipKmerCount = 0;
		this.isExact = pair.isBreakendExact();
		this.be = pair.getBreakendSummary();
		this.category = pair.getEvidenceSource().getSourceCategory();
		this.kmers = new PackedKmerList(k, bases, quals, shouldReverseComplement, shouldReverseComplement);
		this.ambiguous = markAmbiguous(k, bases);
	}
	public VariantEvidence(int k, SingleReadEvidence softClipEvidence, LinearGenomicCoordinate lgc) {
		SAMRecord read = softClipEvidence.getSAMRecord();
		this.evidenceID = softClipEvidence.getEvidenceID();
		boolean shouldReverseComplement = false;		
		byte[] bases = read.getReadBases();
		byte[] quals = read.getBaseQualities();
		
		// calculate anchor info 
		List<CigarElement> elements = read.getCigar().getCigarElements();
		int startClipLength = SAMRecordUtil.getStartSoftClipLength(elements);
		int endClipLength = SAMRecordUtil.getEndSoftClipLength(elements);
		int readLength = bases.length;
		// SSSSMMMMMMMSSSSSSS len=18, start=4, ref=7, end=7, k=4
		// SSSS      |      | ignored kmer
		//  SSSM     |      | ignored kmer
		//   SSMM    |      | ignored kmer
		//    SMMM   |      | ignored kmer
		//	   MMMM  |      | ref kmer
		//      MMMM |      | ref kmer
		//       MMMM|      | ref kmer
		//        MMMM      | ref kmer     <- anchor kmer
		//         MMMS     | branch kmer
		//          MMSS    | branch kmer
		//           MMSS   | branch kmer
		//            MSSS  | branch kmer
		//            |SSSS | branch kmer
		//            | SSSS| branch kmer
		//            |  SSSS branch kmer
		//            ^
		//         anchor pos
		this.referenceKmerStartOffset = startClipLength;
		this.referenceKmerEndOffset = readLength - k - endClipLength + 1;
		int chrPos;
		if (softClipEvidence.getBreakendSummary().direction == BreakendDirection.Forward) {
			this.startSkipKmerCount = startClipLength;
			this.endSkipKmerCount = 0;
			chrPos = read.getUnclippedEnd() - readLength + 1;
		} else {
			this.startSkipKmerCount = 0;
			this.endSkipKmerCount = endClipLength;
			chrPos = read.getUnclippedStart();
		}
		this.expectedAnchorPos = lgc.getLinearCoordinate(softClipEvidence.getBreakendSummary().referenceIndex, chrPos);
		this.isExact = softClipEvidence.isBreakendExact();
		this.be = softClipEvidence.getBreakendSummary();
		this.category = softClipEvidence.getEvidenceSource().getSourceCategory();
		this.kmers = new PackedKmerList(k, bases, quals, shouldReverseComplement, shouldReverseComplement);
		this.ambiguous = markAmbiguous(k, bases);
	}
	private static BitSet markAmbiguous(int k, byte[] bases) {
		BitSet lookup = new BitSet(Math.max(0,bases.length - k + 1));
		for (int i = 0; i < bases.length; i++) {
			if (!SequenceUtil.isValidBase(bases[i])) {
				lookup.set(Math.max(0, i - k + 1), Math.min(lookup.size(), i + 1));
			}
		}
		return lookup;
	}
	public PackedKmerList getKmers() {
		return kmers;
	}
	public int getReferenceKmerCount() {
		return referenceKmerEndOffset - referenceKmerStartOffset;
	}
	public boolean isSkippedKmer(int readKmerOffset) {
		return readKmerOffset < startSkipKmerCount || readKmerOffset >= kmerCount() - endSkipKmerCount; 
	}
	/**
	 * Determines whether this kmer fully supports the reference
	 * @param readKmerOffset kmer index of read
	 * @return true first kmer of this read to be included in the graph, false otherwise 
	 */
	public boolean isReferenceKmer(int readKmerOffset) {
		return readKmerOffset >= referenceKmerStartOffset && readKmerOffset < referenceKmerEndOffset;
	}
	/**
	 * Linear genomic anchor position of reference of the first base of this kmer
	 * as if the evidence supported the reference allele 
	 * @param readKmerOffset kmer index of read
	 * @return inferred reference anchor of this kmer
	 */
	public long getExpectedLinearPosition(int readKmerOffset) {
		return expectedAnchorPos + readKmerOffset;
	}
	private int kmerCount() { return kmers.length(); }
	/**
	 * @return true if this evidence is directly anchored to the reference (ie: direct breakend evidence),
	 * false if anchored by the mate 
	 */
	public boolean isDirectlyAnchoredToReference() { return isExact; }
	public BreakendSummary getBreakend() { return be; }
	public int getCategory() { return category; }
	/**
	 * Gets the underlying evidence
	 * @return evidence
	 */
	public String getEvidenceID() {
		return evidenceID;
	}
	public boolean containsAmbiguousBases(int readKmerOffset) {
		return ambiguous.get(readKmerOffset);
	}
}
