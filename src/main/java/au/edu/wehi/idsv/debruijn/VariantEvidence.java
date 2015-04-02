package au.edu.wehi.idsv.debruijn;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

import java.util.List;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.LinearGenomicCoordinate;
import au.edu.wehi.idsv.NonReferenceReadPair;
import au.edu.wehi.idsv.SoftClipEvidence;
import au.edu.wehi.idsv.sam.SAMRecordUtil;

/**
 * Helper wrapper exposing information both the raw evidence read and
 * kmer-specific information for that read.
 *   
 * @author Daniel Cameron
 *
 */
public class VariantEvidence {
	private final boolean shouldReverseComplement;
	private final byte[] bases;
	private final byte[] quals;
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
	private final int k;
	private final DirectedEvidence evidence;
	
	public VariantEvidence(int k, NonReferenceReadPair pair, LinearGenomicCoordinate lgc) {
		this.k = k;
		this.evidence = pair;
		this.shouldReverseComplement = !pair.onExpectedStrand();		
		this.bases = pair.getNonReferenceRead().getReadBases();
		this.quals = pair.getNonReferenceRead().getBaseQualities();
		int chrPos;
		if (evidence.getBreakendSummary().direction == BreakendDirection.Forward) {
			chrPos = pair.getNonReferenceRead().getUnclippedStart() + pair.getEvidenceSource().getExpectedFragmentSize() - bases.length;
		} else {
			chrPos = pair.getNonReferenceRead().getUnclippedEnd() - pair.getEvidenceSource().getExpectedFragmentSize() + 1;
		}
		this.expectedAnchorPos = lgc.getLinearCoordinate(evidence.getBreakendSummary().referenceIndex, chrPos);
		this.referenceKmerStartOffset = -1;
		this.referenceKmerEndOffset = -1;
		this.startSkipKmerCount = 0;
		this.endSkipKmerCount = 0;
	}
	public VariantEvidence(int k, SoftClipEvidence softClipEvidence, LinearGenomicCoordinate lgc) {
		SAMRecord read = softClipEvidence.getSAMRecord();
		this.k = k;
		this.evidence = softClipEvidence;
		this.shouldReverseComplement = false;		
		this.bases = read.getReadBases();
		this.quals = read.getBaseQualities();
		
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
		this.referenceKmerEndOffset = readLength - k - endClipLength;
		int chrPos;
		if (evidence.getBreakendSummary().direction == BreakendDirection.Forward) {
			this.startSkipKmerCount = startClipLength;
			this.endSkipKmerCount = 0;
			chrPos = read.getUnclippedEnd() - readLength + 1;
		} else {
			this.startSkipKmerCount = 0;
			this.endSkipKmerCount = endClipLength;
			chrPos = read.getUnclippedStart();
		}
		this.expectedAnchorPos = lgc.getLinearCoordinate(evidence.getBreakendSummary().referenceIndex, chrPos);
	}
	public ReadKmerIterable getReadKmers() {
		return new ReadKmerIterable(k, bases, quals, shouldReverseComplement, shouldReverseComplement);
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
	public int kmerCount() { return bases.length - k + 1; }
	/**
	 * @return true if this evidence is directly anchored to the reference (ie: direct breakend evidence),
	 * false if anchored by the mate 
	 */
	public boolean isDirectlyAnchoredToReference() { return evidence.isBreakendExact(); }
	/**
	 * Gets the underlying evidence
	 * @return evidence
	 */
	public DirectedEvidence getDirectedEvidence() {
		return evidence;
	}
}
