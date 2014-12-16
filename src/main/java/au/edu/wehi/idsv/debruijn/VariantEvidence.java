package au.edu.wehi.idsv.debruijn;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

import java.util.List;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.NonReferenceReadPair;
import au.edu.wehi.idsv.SoftClipEvidence;
import au.edu.wehi.idsv.sam.SAMRecordUtil;

import com.google.common.collect.Lists;

/**
 * Helper wrapper exposing information both the raw evidence read and
 * kmer-specific information for that read.
 *   
 * @author Daniel Cameron
 *
 */
public class VariantEvidence {
	private static boolean shouldReverse(SAMRecord record, boolean isAnchored, BreakendDirection graphDirection) {
		boolean reverseDirection = graphDirection == BreakendDirection.Backward;
		return reverseDirection ^ isReverseComplimentRequiredForPositiveStrandReadout(record, isAnchored);
	}
	private static boolean isReverseComplimentRequiredForPositiveStrandReadout(SAMRecord record, boolean isAnchored) {
		if (isAnchored && record.getReadUnmappedFlag()) throw new RuntimeException("Sanity check failure: anchored read is unmapped");
		if (isAnchored) return false;
		if (!record.getReadPairedFlag()) throw new RuntimeException("Sanity check failure: Cannot determine expected orientation of unanchored unpaired read");
		if (record.getMateUnmappedFlag()) throw new RuntimeException("Sanity check failure: Cannot determine expected orientation of unanchored read with unmapped mate");
		// We expect FR pairing orientation
		// Mate anchored on negative strand = need positive strand mate
		// 	-> reverse if mate+ & read+ or mate- & read- === mate XOR read
		return record.getReadNegativeStrandFlag() == record.getMateNegativeStrandFlag();
	}
	private static boolean shouldCompliment(SAMRecord record, boolean anchored) {
		return isReverseComplimentRequiredForPositiveStrandReadout(record, anchored);
	}
	private final boolean isReversed;
	private final boolean isComplemented;
	private final SAMRecord record;
	private final int startSkipKmerCount;
	private final int referenceKmerCount;
	private final int referenceKmerAnchorPosition;
	private final int mateAnchorPosition;
	private final int k;
	private final BreakendDirection direction;
	private final DirectedEvidence evidence;
	
	public static VariantEvidence createRemoteReadEvidence(BreakendDirection graphDirection, int k, NonReferenceReadPair pair) {
		if (graphDirection != pair.getBreakendSummary().direction) throw new RuntimeException("Sanity check failure: local read pair evidence direction does not match de bruijn graph direction");
		int mateAnchorPosition = graphDirection == BreakendDirection.Forward ? pair.getBreakendSummary().start: pair.getBreakendSummary().end;
		return new VariantEvidence(pair, graphDirection, k, pair.getNonReferenceRead(), mateAnchorPosition);
	}
	public static VariantEvidence createSoftClipEvidence(BreakendDirection graphDirection, int k, SoftClipEvidence read) {
		if (graphDirection != read.getBreakendSummary().direction) throw new RuntimeException("Sanity check failure: soft clip direction does not match de bruijn graph direction");
		return new VariantEvidence(read, graphDirection, k, read.getSAMRecord());
	}
	/**
	 * Creates unanchored De Bruijn graph read evidence
	 */
	private VariantEvidence(DirectedEvidence evidence, BreakendDirection graphDirection, int k, SAMRecord record, int mateAnchorPosition) {
		this.evidence = evidence;
		this.direction = graphDirection;
		this.k = k;
		this.record = record;
		this.isReversed = shouldReverse(record, false, graphDirection);
		this.isComplemented = shouldCompliment(record, false);
		this.startSkipKmerCount = 0; // ignore the initial kmers that are soft clipped in the incorrect direction
		this.referenceKmerCount = 0;
		this.referenceKmerAnchorPosition = Integer.MIN_VALUE;
		this.mateAnchorPosition = mateAnchorPosition;
	}
	/**
	 * Creates anchored De Bruijn graph read evidence
	 */
	private VariantEvidence(DirectedEvidence evidence, BreakendDirection graphDirection, int k, SAMRecord record) {
		this.evidence = evidence;
		this.direction = graphDirection;
		this.k = k;
		this.record = record;
		this.isReversed = shouldReverse(record, true, graphDirection);
		this.isComplemented = shouldCompliment(record, true);
		List<CigarElement> elements = isReversed ? Lists.reverse(record.getCigar().getCigarElements()) : record.getCigar().getCigarElements();
		int startClipLength = SAMRecordUtil.getStartSoftClipLength(elements);
		int endClipLength = SAMRecordUtil.getEndSoftClipLength(elements);
		int readLength = record.getReadLength();
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
		this.startSkipKmerCount = startClipLength; // ignore the initial kmers that are soft clipped in the incorrect direction
		int branchKmers = endClipLength;
		int totalKmers = readLength - k + 1;
		this.referenceKmerCount = totalKmers - startSkipKmerCount - branchKmers;
		this.referenceKmerAnchorPosition = isReversed ? record.getAlignmentStart() : record.getAlignmentEnd();
		this.mateAnchorPosition = Integer.MIN_VALUE;
	}
	public boolean isReversed() {
		return isReversed;
	}
	public boolean isComplemented() {
		return isComplemented;
	}
	public int getStartSkipKmerCount() {
		return startSkipKmerCount;
	}
	public int getReferenceKmerCount() {
		return referenceKmerCount;
	}
	/**
	 * Genomic position of base of reference-supporting kmer closest to the putative structural variation 
	 * @return
	 */
	public int getReferenceKmerAnchorPosition() {
		if (referenceKmerAnchorPosition == Integer.MIN_VALUE) throw new IllegalStateException(String.format("Not anchored evidence."));
		return referenceKmerAnchorPosition;
	}
	/**
	 * Genomic position of mapped mate read
	 * 
	 * WARNING: assumes FR orientation
	 * @return
	 */
	public int getMateAnchorPosition() {
		if (mateAnchorPosition == Integer.MIN_VALUE) throw new IllegalStateException(String.format("Not unanchored evidence."));
		return mateAnchorPosition;
	}
	public SAMRecord getSAMRecord() {
		return record;
	}
	public BreakendDirection getDirection() {
		return direction;
	}
	/**
	 * Determines whether this is the first kmer of this read to be included in the graph
	 * @param readKmerOffset kmer index of read
	 * @return true first kmer of this read to be included in the graph, false otherwise 
	 */
	public boolean isFirstKmer(int readKmerOffset) {
		return readKmerOffset == startSkipKmerCount;
	}
	public boolean isSkippedKmer(int readKmerOffset) {
		return readKmerOffset >= firstSkippedKmerOffset() && readKmerOffset <= lastSkippedKmerOffset();
	}
	/**
	 * Determines whether this kmer fully supports the reference
	 * @param readKmerOffset kmer index of read
	 * @return true first kmer of this read to be included in the graph, false otherwise 
	 */
	public boolean isReferenceKmer(int readKmerOffset) {
		return readKmerOffset >= firstReferenceKmerOffset() && readKmerOffset <= lastReferenceKmerOffset();
	}
	/**
	 * Determines whether this kmer provides evidence for a variant
	 * @param readKmerOffset kmer index of read
	 * @return
	 */
	public boolean isVariantKmer(int readKmerOffset) {
		return readKmerOffset >= firstVariantKmerOffset();
	}
	/**
	 * How many bases of this kmer support the reference
	 * @param readKmerOffset kmer index of read
	 * @return number of kmer bases supporting the reference
	 */
	public int basesSupportingReference(int readKmerOffset) {
		if (referenceKmerCount == 0) return 0;
		if (isReferenceKmer(readKmerOffset)) return k;
		if (isVariantKmer(readKmerOffset)) return Math.max(0, startSkipKmerCount + referenceKmerCount + k - readKmerOffset - 1);
		// isSkippedKmer == true
		return Math.max(0, readKmerOffset + k - startSkipKmerCount);
	}
	public boolean isReferenceAnchor(int readKmerOffset) {
		return referenceKmerCount > 0 && readKmerOffset == lastReferenceKmerOffset();
	}
	public boolean isMateAnchor(int readKmerOffset) {
		return referenceKmerCount == 0 && isFirstKmer(readKmerOffset);
	}
	/**
	 * Genomic anchor position of reference kmer assuming read is directly aligned
	 * to the reference between the anchor and this position
	 * @param readKmerOffset kmer index of read
	 * @return inferred reference anchor of this kmer
	 */
	public int getInferredReferencePosition(int readKmerOffset) {
		if (referenceKmerAnchorPosition == Integer.MIN_VALUE) throw new IllegalStateException(String.format("Not anchored evidence."));
		int offset = lastReferenceKmerOffset() - readKmerOffset;
		return referenceKmerAnchorPosition - offset * (direction == BreakendDirection.Forward ? 1 : -1);
	}
	/**
	 * Returns the inferred genomic position of the first base of the kmer.
	 * @param readKmerOffset kmer index of read
	 * @return genomic position of first base
	 */
	public int getReferenceStartingPosition(int readKmerOffset) {
		if (basesSupportingReference(readKmerOffset) == 0 && isVariantKmer(readKmerOffset)) throw new IllegalArgumentException("no reference bases in kmer");
		int offsetFromAnchor = readKmerOffset - lastReferenceKmerOffset();
		return referenceKmerAnchorPosition - (offsetFromAnchor - k + 1) * genomicPositionDeltaPerKmer();
	}
	/**
	 * Change in genomic position when moving to the next kmer
	 * @return 1 if we move to the next genomic position for the next kmer, -1 if we go backwards over the chromosome
	 */
	private int genomicPositionDeltaPerKmer() {
		if (referenceKmerCount > 0) return isReversed ? 1 : -1;
		throw new RuntimeException("Unable to determine genomic position delta for unanchored reads");
	}
	public int firstSkippedKmerOffset() { return startSkipKmerCount == 0 ? -1 : 0; }
	public int lastSkippedKmerOffset() { return startSkipKmerCount - 1; }
	public int firstReferenceKmerOffset() { return referenceKmerCount == 0 ? -1 : startSkipKmerCount; }
	public int lastReferenceKmerOffset() { return referenceKmerCount == 0 ? -1 : (startSkipKmerCount + referenceKmerCount - 1); }
	public int firstVariantKmerOffset() { return startSkipKmerCount + referenceKmerCount; }
	public int lastVariantKmerOffset() { return lastKmerOffset(); }
	public int lastKmerOffset() { return kmerCount() - 1; }
	public int kmerCount() { return record.getReadLength() - k + 1; }
	/**
	 * @return true if this evidence is directly anchored to the reference (ie: direct breakend evidence),
	 * false if anchored by the mate 
	 */
	public boolean isDirectlyAnchoredToReference() { return referenceKmerAnchorPosition != Integer.MIN_VALUE; }
	/**
	 * Gets the underlying evidence
	 * @return evidence
	 */
	public DirectedEvidence getDirectedEvidence() {
		return evidence;
	}
}
