/**
 * 
 */
package au.edu.wehi.idsv.sam;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamPairUtil;
import htsjdk.samtools.SamPairUtil.PairOrientation;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.SequenceUtil;

import java.util.Arrays;
import java.util.List;

/**
 * @author Daniel Cameron
 *
 */
public class SAMRecordUtil {
	/**
	 * Determines whether the given record is soft-clipped
	 * @param aln
	 * @return true if soft-clipped, false otherwise
	 */
	public static boolean isAlignmentSoftClipped(SAMRecord aln) {
		return isSoftClipLengthAtLeast(aln, 1);
	}
	/**
	 * Determines whether the given record has a soft clip at least the given length
	 * @param aln
	 * @param length minimum length of softclip
	 * @return true if soft-clipped and at least length
	 */
	public static boolean isSoftClipLengthAtLeast(SAMRecord aln, int length) {
		return !aln.getReadUnmappedFlag() && (getStartSoftClipLength(aln) >= length || getEndSoftClipLength(aln) >= length);
	}
	public static int getStartSoftClipLength(SAMRecord aln) {
		Cigar cigar = aln.getCigar();
		if (cigar == null) return 0;
		return getStartSoftClipLength(cigar.getCigarElements());
	}
	public static int getStartSoftClipLength(List<CigarElement> elements) {
		if (elements == null) return 0;
		int i = 0;
		while (i < elements.size() && (elements.get(i).getOperator() == CigarOperator.HARD_CLIP || elements.get(i).getOperator() == CigarOperator.SOFT_CLIP)) {
			if (elements.get(i).getOperator() == CigarOperator.SOFT_CLIP) return elements.get(i).getLength();
			i++;
		}
		return 0;
	}
	public static int getEndSoftClipLength(SAMRecord aln) {
		Cigar cigar = aln.getCigar();
		if (cigar == null) return 0;
		return getEndSoftClipLength(cigar.getCigarElements());
	}
	public static int getEndSoftClipLength(List<CigarElement> elements) {
		if (elements == null) return 0;
		int i = elements.size() - 1;
		while (i > 0 && (elements.get(i).getOperator() == CigarOperator.HARD_CLIP || elements.get(i).getOperator() == CigarOperator.SOFT_CLIP)) {
			if (elements.get(i).getOperator() == CigarOperator.SOFT_CLIP) return elements.get(i).getLength();
			i--;
		}
		return 0;
	}
	public static byte[] getStartSoftClipBases(SAMRecord record) {
		byte[] seq = record.getReadBases();
		if (seq == null) return null;
		if (seq == SAMRecord.NULL_SEQUENCE) return null;
		seq = Arrays.copyOfRange(seq, 0, getStartSoftClipLength(record));
		return seq;
	}
	public static byte[] getEndSoftClipBases(SAMRecord record) {
		byte[] seq = record.getReadBases();
		if (seq == null) return null;
		if (seq == SAMRecord.NULL_SEQUENCE) return null;
		seq = Arrays.copyOfRange(seq, record.getReadLength() - getEndSoftClipLength(record), record.getReadLength());
		return seq;
	}
	public static byte[] getStartSoftClipBaseQualities(SAMRecord record) {
		byte[] seq = record.getBaseQualities();
		if (seq == null) return null;
		if (seq == SAMRecord.NULL_QUALS) return null;
		seq = Arrays.copyOfRange(seq, 0, getStartSoftClipLength(record));
		return seq;
	}
	public static byte[] getEndSoftClipBaseQualities(SAMRecord record) {
		byte[] seq = record.getBaseQualities();
		if (seq == null) return null;
		if (seq == SAMRecord.NULL_QUALS) return null;
		seq = Arrays.copyOfRange(seq, record.getReadLength() - getEndSoftClipLength(record), record.getReadLength());
		return seq;
	}
	/**
	 * Determines whether this read is part of a non-reference read pair
	 * @param record SAMRecord to check
	 * @return true if part of non-reference pair, false otherwise
	 */
	public static boolean isPartOfNonReferenceReadPair(SAMRecord record) {
		return record != null &&
				record.getReadPairedFlag() &&
				!record.getProperPairFlag() &&
				// at least one of the reads should be mapped
				!(record.getReadUnmappedFlag() && record.getMateUnmappedFlag());
	}
	public static boolean isDiscordantPairMember(SAMRecord record) {
		return record != null &&
				record.getReadPairedFlag() &&
				!record.getProperPairFlag() &&
				!record.getReadUnmappedFlag() &&
				!record.getMateUnmappedFlag();
	}
	public static boolean isAnchoredPairMember(SAMRecord record) {
		return record != null &&
				record.getReadPairedFlag() &&
				// only one read in the pair is mapped
				(record.getReadUnmappedFlag() ^ record.getMateUnmappedFlag());
	}
	public static void ensureNmTag(ReferenceSequenceFileWalker refFileWalker, SAMRecord record) {
		if (record == null) return;
		if (record.getReadBases() == null) return;
		if (record.getReadBases() == SAMRecord.NULL_SEQUENCE) return;
		if (record.getIntegerAttribute(SAMTag.NM.name()) != null) return;
		if (record.getReadUnmappedFlag()) return;
        final ReferenceSequence refSequence = refFileWalker.get(record.getReferenceIndex());
        final int actualNucleotideDiffs = SequenceUtil.calculateSamNmTag(record, refSequence.getBases(), 0, false);
        record.setAttribute(SAMTag.NM.name(), actualNucleotideDiffs);
	}
	public static void ensureNmTag(ReferenceSequenceFile ref, SAMRecord record) {
		if (record == null) return;
		if (record.getReadBases() == null) return;
		if (record.getReadBases() == SAMRecord.NULL_SEQUENCE) return;
		if (record.getIntegerAttribute(SAMTag.NM.name()) != null) return;
		if (record.getReadUnmappedFlag()) return;
        final ReferenceSequence refSequence = ref.getSequence(ref.getSequenceDictionary().getSequence(record.getReferenceIndex()).getSequenceName());
        final int actualNucleotideDiffs = SequenceUtil.calculateSamNmTag(record, refSequence.getBases(), 0, false);
        record.setAttribute(SAMTag.NM.name(), actualNucleotideDiffs);
	}
	public static int getMaxReferenceBaseQual(SAMRecord r) {
		byte[] qual = r.getBaseQualities();
		if (qual == null || qual == SAMRecord.NULL_QUALS) return 0;
		int max = 0;
		for (int i = getStartSoftClipLength(r); i < qual.length - getEndSoftClipLength(r); i++) {
			max = Math.max(max, qual[i]);
		}
		return max;
	}
	public static int getTotalReferenceBaseQual(SAMRecord r) {
		byte[] qual = r.getBaseQualities();
		if (qual == null || qual == SAMRecord.NULL_QUALS) return 0;
		int max = 0;
		for (int i = getStartSoftClipLength(r); i < qual.length - getEndSoftClipLength(r); i++) {
			max += qual[i];
		}
		return max;
	}
	/**
	 * Determines whether the given read pair mapping could be a valid fragment
	 * @param record paired read
	 * @param expectedOrientation expected orientation
	 * @return true if the given read pair could come from a reference-supporting fragment 
	 */
	public static boolean couldBeProperPair(SAMRecord record, PairOrientation expectedOrientation) {
		return record.getReadPairedFlag()
				&& !record.getReadUnmappedFlag()
				&& !record.getMateUnmappedFlag()
				&& record.getReferenceIndex() == record.getMateReferenceIndex()
				&& SamPairUtil.getPairOrientation(record) == expectedOrientation;
	}	
}
