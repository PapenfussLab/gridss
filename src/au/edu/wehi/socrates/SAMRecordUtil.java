/**
 * 
 */
package au.edu.wehi.socrates;

import java.util.List;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFileWalker;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMTag;
import net.sf.samtools.util.SequenceUtil;

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
		return !aln.getReadUnmappedFlag() && (getStartSoftClipLength(aln) > 0 || getEndSoftClipLength(aln) > 0);
	}
	public static int getStartSoftClipLength(SAMRecord aln) {
		Cigar cigar = aln.getCigar();
		if (cigar == null) return 0;
		List<CigarElement> elements = cigar.getCigarElements();
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
		List<CigarElement> elements = cigar.getCigarElements();
		if (elements == null) return 0;
		int i = elements.size() - 1;
		while (i > 0 && (elements.get(i).getOperator() == CigarOperator.HARD_CLIP || elements.get(i).getOperator() == CigarOperator.SOFT_CLIP)) {
			if (elements.get(i).getOperator() == CigarOperator.SOFT_CLIP) return elements.get(i).getLength();
			i--;
		}
		return 0;
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
        final ReferenceSequence refSequence = refFileWalker.get(record.getReferenceIndex());
        final int actualNucleotideDiffs = SequenceUtil.calculateSamNmTag(record, refSequence.getBases(), 0, false);
        record.setAttribute(SAMTag.NM.name(), actualNucleotideDiffs);
	}
}
