/**
 * 
 */
package au.edu.wehi.idsv.sam;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamPairUtil.PairOrientation;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SequenceUtil;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.Pair;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.alignment.AlignerFactory;
import au.edu.wehi.idsv.alignment.Alignment;
import au.edu.wehi.idsv.picard.ReferenceLookup;

import com.google.common.collect.Lists;

/**
 * @author Daniel Cameron
 *
 */
public class SAMRecordUtil {
	private static final Log log = Log.getInstance(SAMRecordUtil.class);
	public static final String FIRST_OF_PAIR_NAME_SUFFIX = "\\1";
	public static final String SECOND_OF_PAIR_NAME_SUFFIX = "\\2";
	/**
	 * Determines whether the given record is soft-clipped
	 * 
	 * @param aln
	 * @return true if soft-clipped, false otherwise
	 */
	public static boolean isAlignmentSoftClipped(SAMRecord aln) {
		return isSoftClipLengthAtLeast(aln, 1);
	}

	/**
	 * Determines whether the given record has a soft clip at least the given
	 * length
	 * 
	 * @param aln
	 * @param length
	 *            minimum length of softclip
	 * @return true if soft-clipped and at least length
	 */
	public static boolean isSoftClipLengthAtLeast(SAMRecord aln, int length) {
		return !aln.getReadUnmappedFlag()
				&& (getStartSoftClipLength(aln) >= length || getEndSoftClipLength(aln) >= length);
	}

	public static int getStartSoftClipLength(SAMRecord aln) {
		Cigar cigar = aln.getCigar();
		if (cigar == null)
			return 0;
		return getStartSoftClipLength(cigar.getCigarElements());
	}

	public static int getStartSoftClipLength(List<CigarElement> elements) {
		if (elements == null)
			return 0;
		int i = 0;
		while (i < elements.size()
				&& (elements.get(i).getOperator() == CigarOperator.HARD_CLIP || elements
						.get(i).getOperator() == CigarOperator.SOFT_CLIP)) {
			if (elements.get(i).getOperator() == CigarOperator.SOFT_CLIP)
				return elements.get(i).getLength();
			i++;
		}
		return 0;
	}

	public static int getEndSoftClipLength(SAMRecord aln) {
		Cigar cigar = aln.getCigar();
		if (cigar == null)
			return 0;
		return getEndSoftClipLength(cigar.getCigarElements());
	}

	public static int getEndSoftClipLength(List<CigarElement> elements) {
		if (elements == null)
			return 0;
		int i = elements.size() - 1;
		while (i > 0
				&& (elements.get(i).getOperator() == CigarOperator.HARD_CLIP || elements
						.get(i).getOperator() == CigarOperator.SOFT_CLIP)) {
			if (elements.get(i).getOperator() == CigarOperator.SOFT_CLIP)
				return elements.get(i).getLength();
			i--;
		}
		return 0;
	}
	
	public static int getSoftClipLength(List<CigarElement> elements, BreakendDirection direction) {
		if (direction == BreakendDirection.Forward) {
			return getEndSoftClipLength(elements);
		} else {
			return getStartSoftClipLength(elements);
		}
	}
	
	public static int getSoftClipLength(SAMRecord aln, BreakendDirection direction) {
		if (direction == BreakendDirection.Forward) {
			return getEndSoftClipLength(aln);
		} else {
			return getStartSoftClipLength(aln);
		}
	}

	public static byte[] getStartSoftClipBases(SAMRecord record) {
		byte[] seq = record.getReadBases();
		if (seq == null)
			return null;
		if (seq == SAMRecord.NULL_SEQUENCE)
			return null;
		seq = Arrays.copyOfRange(seq, 0, getStartSoftClipLength(record));
		return seq;
	}

	public static byte[] getEndSoftClipBases(SAMRecord record) {
		byte[] seq = record.getReadBases();
		if (seq == null)
			return null;
		if (seq == SAMRecord.NULL_SEQUENCE)
			return null;
		seq = Arrays.copyOfRange(seq, record.getReadLength()
				- getEndSoftClipLength(record), record.getReadLength());
		return seq;
	}

	public static byte[] getStartSoftClipBaseQualities(SAMRecord record) {
		byte[] seq = record.getBaseQualities();
		if (seq == null)
			return null;
		if (seq == SAMRecord.NULL_QUALS)
			return null;
		seq = Arrays.copyOfRange(seq, 0, getStartSoftClipLength(record));
		return seq;
	}

	public static byte[] getEndSoftClipBaseQualities(SAMRecord record) {
		byte[] seq = record.getBaseQualities();
		if (seq == null)
			return null;
		if (seq == SAMRecord.NULL_QUALS)
			return null;
		seq = Arrays.copyOfRange(seq, record.getReadLength()
				- getEndSoftClipLength(record), record.getReadLength());
		return seq;
	}

	public static void ensureNmTag(ReferenceSequenceFileWalker refFileWalker,
			SAMRecord record) {
		if (record == null)
			return;
		if (record.getReadBases() == null)
			return;
		if (record.getReadBases() == SAMRecord.NULL_SEQUENCE)
			return;
		if (record.getIntegerAttribute(SAMTag.NM.name()) != null)
			return;
		if (record.getReadUnmappedFlag())
			return;
		final ReferenceSequence refSequence = refFileWalker.get(record
				.getReferenceIndex());
		final int actualNucleotideDiffs = SequenceUtil.calculateSamNmTag(
				record, refSequence.getBases(), 0, false);
		record.setAttribute(SAMTag.NM.name(), actualNucleotideDiffs);
	}

	public static void ensureNmTag(ReferenceSequenceFile ref, SAMRecord record) {
		if (record == null)
			return;
		if (record.getReadBases() == null)
			return;
		if (record.getReadBases() == SAMRecord.NULL_SEQUENCE)
			return;
		if (record.getIntegerAttribute(SAMTag.NM.name()) != null)
			return;
		if (record.getReadUnmappedFlag())
			return;
		final ReferenceSequence refSequence = ref.getSequence(ref
				.getSequenceDictionary()
				.getSequence(record.getReferenceIndex()).getSequenceName());
		final int actualNucleotideDiffs = SequenceUtil.calculateSamNmTag(
				record, refSequence.getBases(), 0, false);
		record.setAttribute(SAMTag.NM.name(), actualNucleotideDiffs);
	}

	public static int getMaxReferenceBaseQual(SAMRecord r) {
		byte[] qual = r.getBaseQualities();
		if (qual == null || qual == SAMRecord.NULL_QUALS)
			return 0;
		int max = 0;
		for (int i = getStartSoftClipLength(r); i < qual.length
				- getEndSoftClipLength(r); i++) {
			max = Math.max(max, qual[i]);
		}
		return max;
	}

	public static int getTotalReferenceBaseQual(SAMRecord r) {
		byte[] qual = r.getBaseQualities();
		if (qual == null || qual == SAMRecord.NULL_QUALS)
			return 0;
		int max = 0;
		for (int i = getStartSoftClipLength(r); i < qual.length
				- getEndSoftClipLength(r); i++) {
			max += qual[i];
		}
		return max;
	}
	/**
	 * Dovetailing reads either either due to an SV
	 * or failure to trim adapters from a fragment smaller than the read length
	 * 
	 *    =======>
	 * <=======
	 * 
	 * @param expectedOrientation
	 *            read pair orientation
	 * @return true if the soft clip is due to a fragment size smaller than the
	 *         read length
	 */
	public static boolean isDovetailing(SAMRecord record1, SAMRecord record2, PairOrientation expectedOrientation, int margin) {
		if (expectedOrientation != PairOrientation.FR) throw new RuntimeException("NYI");
		return !record1.getReadUnmappedFlag()
				&& !record2.getReadUnmappedFlag()
				&& record1.getReferenceIndex() == record2.getReferenceIndex()
				&& record1.getReadNegativeStrandFlag() != record2.getReadNegativeStrandFlag() // FR
				&& overlap(record1, record2)
				&& Math.abs(record1.getAlignmentStart()
						- record2.getAlignmentStart()) <= margin
				&& Math.abs(record1.getAlignmentEnd()
						- record2.getAlignmentEnd()) <= margin;
	}
	public static boolean overlap(SAMRecord r1, SAMRecord r2) {
		boolean result = r1 != null
				&& r2 != null
				&& !r1.getReadUnmappedFlag()
				&& !r2.getReadUnmappedFlag()
				&& r1.getReferenceIndex() == r2.getReferenceIndex()
				&& ((r1.getAlignmentStart() >= r2.getAlignmentStart() && r1
						.getAlignmentStart() <= r2.getAlignmentEnd()) || (r2
						.getAlignmentStart() >= r1.getAlignmentStart() && r2
						.getAlignmentStart() <= r1.getAlignmentEnd())); 
		return result;
	}
	/**
	 * Estimates the size of sequenced fragment
	 * @param record
	 * @param expectedOrientation 
	 * @return
	 */
	public static int estimateFragmentSize(SAMRecord record, PairOrientation expectedOrientation) {
		if (expectedOrientation != PairOrientation.FR) throw new RuntimeException("NYI");
		if (record.getReadUnmappedFlag() ||
				record.getReadUnmappedFlag() ||
				record.getMateUnmappedFlag() ||
				record.getReferenceIndex() != record.getMateReferenceIndex() ||
				// FR assumption
				record.getReadNegativeStrandFlag() == record.getMateNegativeStrandFlag()) {
			return 0;
		}
		// Assuming FR orientation, adapter sequences have been removed, and no SCs
		if (record.getReadNegativeStrandFlag()) {
			// <--record
			return record.getUnclippedEnd() - record.getMateAlignmentStart() + 1;
		} else {
			// record--> 
			// this assumes the mate is fully mapped to the reference (best guess we can make without mate cigar)
			return record.getMateAlignmentStart() + record.getReadLength() - record.getUnclippedStart();
		}
	}
	/**
	 * Conservative estimate as to whether the reads overlap due to small fragment size
	 * @param local
	 * @param remote
	 * @return number possible breakpoints between the read pair mapped in the expected orientation,
	 *  or Integer.MAX_VALUE if placement is not as expected
	 */
	public static boolean estimatedReadsOverlap(SAMRecord record, PairOrientation expectedOrientation, int minExpectedAnchorBases) {
		if (expectedOrientation != PairOrientation.FR) throw new RuntimeException("NYI");
		if (record.getReadUnmappedFlag() ||
				record.getReadUnmappedFlag() ||
				record.getMateUnmappedFlag() ||
				record.getReferenceIndex() != record.getMateReferenceIndex() ||
				// FR assumption
				record.getReadNegativeStrandFlag() == record.getMateNegativeStrandFlag()) {
			return false;
		}
		// Assuming FR orientation, adapter sequences have been removed, and no SCs
		if (record.getReadNegativeStrandFlag()) {
			//   <-----	record
			// ----->	mate
			if (record.getMateAlignmentStart() > record.getAlignmentStart()) return false; // not FR
			// we don't know where the mate actually ends
			// we can't assume = read length as that will
			// remove our pair when it spans a deletion
			// and our mate has an inner soft clip or non-reference CIGAR
			// so we make a conservative guess based on the min bases aligners
			// will return a alignment for
			int offset = Math.min(minExpectedAnchorBases, record.getReadLength()) - 1;
			return record.getMateAlignmentStart() + offset >= record.getAlignmentStart() &&
					record.getMateAlignmentStart() + offset <= record.getAlignmentEnd();
			
		} else {
			// record ----->
			//           <-----
			return record.getMateAlignmentStart() >= record.getAlignmentStart() &&
					record.getMateAlignmentStart() <= record.getAlignmentEnd();
		}
	}
	/**
	 * Estimates the size of sequenced fragment
	 * @param record
	 * @return
	 */
	public static int calculateFragmentSize(SAMRecord record1, SAMRecord record2, PairOrientation expectedOrientation) {
		if (expectedOrientation != PairOrientation.FR) throw new RuntimeException("NYI");
		if (record1.getReadUnmappedFlag() ||
			record2.getReadUnmappedFlag() ||
			record1.getReferenceIndex() != record2.getReferenceIndex() ||
			// FR assumption
			record1.getReadNegativeStrandFlag() == record2.getReadNegativeStrandFlag()) {
			return 0;
		}
		// Assuming FR orientation and adapter sequences have been removed 
		if (record1.getReadNegativeStrandFlag()) {
			// <--record
			return record1.getUnclippedEnd() - record2.getUnclippedStart() + 1;
		} else {
			return record2.getUnclippedEnd() - record1.getUnclippedStart() + 1;
		}
	}
	/**
	 * Updates the two reads to indicate the form a read pair
	 * @param first first in pair
	 * @param second second in pair
	 */
	public static void pairReads(SAMRecord first, SAMRecord second) {
		// SAMv1 S2.4
		if (first.getReadUnmappedFlag() && !second.getReadUnmappedFlag()) {
			first.setReferenceIndex(second.getReferenceIndex());
			first.setAlignmentStart(second.getAlignmentStart());
		} else if (second.getReadUnmappedFlag() && !first.getReadUnmappedFlag()) {
			second.setReferenceIndex(first.getReferenceIndex());
			second.setAlignmentStart(first.getAlignmentStart());
		}
		first.setReadPairedFlag(true);
		first.setFirstOfPairFlag(true);
		first.setSecondOfPairFlag(false);
		second.setReadPairedFlag(true);
		second.setFirstOfPairFlag(false);
		second.setSecondOfPairFlag(true);
		first.setMateUnmappedFlag(second.getReadUnmappedFlag());
		first.setMateReferenceIndex(second.getReferenceIndex());
		first.setMateAlignmentStart(second.getAlignmentStart());
		first.setMateNegativeStrandFlag(second.getReadNegativeStrandFlag());
		second.setMateUnmappedFlag(first.getReadUnmappedFlag());
		second.setMateReferenceIndex(first.getReferenceIndex());
		second.setMateAlignmentStart(first.getAlignmentStart());
		second.setMateNegativeStrandFlag(first.getReadNegativeStrandFlag());
		second.setMateUnmappedFlag(first.getReadUnmappedFlag());
		second.setReadName(first.getReadName());
	}
	/**
	 * Temporary helper clone that doesn't throw CloneNotSupportedException
	 * @param r record to clone
	 * @return clone of SAMRecord
	 */
	public static SAMRecord clone(SAMRecord r) {
		try {
			return (SAMRecord)r.clone();
		} catch (CloneNotSupportedException e) {
			throw new SAMException(e);
		}
	}
	/**
	 * Determine whether the two reads are the same. This is useful
	 * for filtering out read pairs which, due to sequencing chemistry
	 * error, both sequence the same end of the fragment.
	 *   
	 * @param r1 read to compare
	 * @param r2 read to compare
	 * @return true if the reads appear to be the same, false otherwise
	 */
	public static boolean areSameRead(SAMRecord r1, SAMRecord r2) {
		// TODO: compare read sequences
		return !r1.getReadUnmappedFlag()
				&& !r2.getReadUnmappedFlag()
				&& r1.getReferenceIndex() == r2.getReferenceIndex()
				&& (r1.getUnclippedStart() == r2.getUnclippedStart() || r1.getUnclippedEnd() == r2.getUnclippedEnd())
				&& r1.getReadNegativeStrandFlag() == r2.getReadNegativeStrandFlag();
	}
	/**
	 * Removes soft-clipped bases from the ends of read
	 */
	public static void trimSoftClips(SAMRecord read, int startCount, int endCount) {
		assert(read.getReadUnmappedFlag() || read.getReadLength() == read.getCigar().getReadLength());
		List<CigarElement> cigar = Lists.newArrayList(read.getCigar().getCigarElements());
		if (startCount > 0) {
			CigarElement sc = cigar.get(0);
			int len = sc.getLength();
			assert(sc.getOperator() == CigarOperator.SOFT_CLIP);
			assert(len >= startCount);
			cigar.remove(0);
			len -= startCount;
			if (len > 0) {
				cigar.add(0, new CigarElement(len, CigarOperator.SOFT_CLIP));
			}
		}
		if (endCount > 0) {
			CigarElement sc = cigar.get(cigar.size() - 1);
			int len = sc.getLength();
			assert(sc.getOperator() == CigarOperator.SOFT_CLIP);
			assert(len >= endCount);
			cigar.remove(cigar.size() - 1);
			len -= endCount;
			if (len > 0) {
				cigar.add(cigar.size(), new CigarElement(len, CigarOperator.SOFT_CLIP));
			}
		}
		int readLength = read.getReadLength();
		read.setCigar(new Cigar(cigar));
		read.setReadBases(Arrays.copyOfRange(read.getReadBases(), startCount, readLength - endCount));
		read.setBaseQualities(Arrays.copyOfRange(read.getBaseQualities(), startCount, readLength - endCount));
		assert(read.getReadLength() == read.getCigar().getReadLength());
	}
	public static void trim(SAMRecord read, int startCount, int endCount) {
		int readLength = read.getReadLength();
		read.setReadBases(Arrays.copyOfRange(read.getReadBases(), startCount, readLength - endCount));
		read.setBaseQualities(Arrays.copyOfRange(read.getBaseQualities(), startCount, readLength - endCount));
		if (read.getCigar() != null && read.getCigar().getCigarElements().size() > 0) {
			read.setAlignmentStart(read.getAlignmentStart() + CigarUtil.offsetOf(read.getCigar(), startCount));
			read.setCigar(CigarUtil.trimReadBases(read.getCigar(), startCount, endCount));
			assert(read.getReadLength() == read.getCigar().getReadLength());
		}
	}
	public static boolean isReferenceAlignment(SAMRecord r) {
		if (r.getCigar() == null) return true;
		return isReferenceAlignment(r.getCigar().getCigarElements());
	}
	public static boolean isReferenceAlignment(List<CigarElement> cigar) {
		if (cigar == null || cigar.size() == 0) return true;
		for (CigarElement e : cigar) {
			if (e.getLength() > 0) {
				switch (e.getOperator()) {
					case M:
					case EQ:
					case X:
					case H:
						break;
					default:
						return false;
				}
			}
		}
		return true;
	}
	/**
	 * Gets the entropy of the read, excluding soft clips
	 * @param read read
	 * @return number of bits of entropy, excluding soft clipped bases
	 */
	public static double alignedEntropy(SAMRecord read) {
		int offset = getStartSoftClipLength(read);
		int length = read.getReadLength() - offset - getEndSoftClipLength(read);
		double bits = au.edu.wehi.idsv.util.SequenceUtil.shannonEntropy(read.getReadBases(), offset, length);
		return bits;
	}
	/**
	 * Gets the entropy of the entire read, regardless of alignment
	 * @param read read
	 * @return number of bits of entropy
	 */
	public static double entropy(SAMRecord read) {
		double bits = au.edu.wehi.idsv.util.SequenceUtil.shannonEntropy(read.getReadBases(), 0, read.getReadLength());
		return bits;
	}
	/**
	 * Splits the read alignment at the given read base position 
	 * @param read read to split
	 * @param readBaseOffset 0-based offset of base to split immediately after
	 * @return left and right split reads 
	 */
	public static Pair<SAMRecord, SAMRecord> splitAfter(SAMRecord read, int readBaseOffset) {
		Pair<Cigar, Cigar> cigars = CigarUtil.splitAfterReadPosition(read.getCigar().getCigarElements(), readBaseOffset);
		SAMRecord left = SAMRecordUtil.clone(read);
		left.setCigar(cigars.getLeft());
		SAMRecord right = SAMRecordUtil.clone(read);
		right.setCigar(cigars.getRight());
		right.setAlignmentStart(read.getAlignmentStart() + CigarUtil.offsetOf(read.getCigar(), CigarUtil.getStartClipLength(cigars.getRight().getCigarElements())));
		return Pair.of(left, right);
	}
	/**
	 * Performs local realignment of the given SAMRecord in a window around the alignment location
	 * @param reference reference genome
	 * @param read read
	 * @param windowSize number of bases to extend window around alignment
	 * @param extendWindowForClippedBases extend window for soft clipped bases
	 * @return copy of realigned read if alignment changed, read if realignment did not change the alignment
	 */
	public static SAMRecord realign(ReferenceLookup reference, SAMRecord read, int windowSize, boolean extendWindowForClippedBases) {
		if (read.getReadUnmappedFlag()) return read;
		SAMSequenceRecord refSeq = reference.getSequenceDictionary().getSequence(read.getReferenceIndex());
		// find reference bounds of read. Negative deletions mean we can't just use
		// getAlignmentStart() and getAlignmentEnd()
		int pos = read.getAlignmentStart();
		int start = pos;
		int end = pos;
		for (CigarElement ce : CigarUtil.decodeNegativeDeletion(read.getCigar().getCigarElements())) {
			if (ce.getOperator().consumesReferenceBases()) {
				pos += ce.getLength();
			}
			start = Math.min(start, pos);
			end = Math.max(end, pos - 1);
		}
		// extend bounds
		start -= windowSize;
		end += windowSize;
		if (extendWindowForClippedBases) {
			start -= getStartSoftClipLength(read);
			end += getEndSoftClipLength(read);
		}
		// don't overrun contig bounds
		start = Math.max(1, start);
		end = Math.min(refSeq.getSequenceLength(), end);
		
		byte[] ass = read.getReadBases();
		byte[] ref = reference.getSubsequenceAt(refSeq.getSequenceName(), start, end).getBases();
		if (ass == null || ref == null || ass.length == 0 || ref.length == 0) {
			return read;
		}
		// defensive checks so we don't crash the JVM if an unexpected character is encountered
		for (int i = 0; i < ass.length; i++) {
			if (!htsjdk.samtools.util.SequenceUtil.isValidBase(ass[i])) {
				ass[i] = 'N';
			}
		}
		for (int i = 0; i < ref.length; i++) {
			if (!htsjdk.samtools.util.SequenceUtil.isValidBase(ref[i])) {
				ref[i] = 'N';
			}
		}
		Alignment alignment =  null;
		try {
			alignment = AlignerFactory.create().align_smith_waterman(ass, ref);        
		} catch (Exception e) {
			// swallow and log alignment error
			log.error(e, String.format("Error aligning %s to %s:%d-%d",
					read.getReadName(),
					refSeq.getSequenceName(), start, end));
		}
		if (alignment == null) {
			return read;
		}
		Cigar cigar = TextCigarCodec.decode(alignment.getCigar());
		
		int alignmentStart = start + alignment.getStartPosition();
		if (alignmentStart == read.getAlignmentStart() && cigar.equals(read.getCigar())) {
			return read;
		}
		SAMRecord copy = SAMRecordUtil.clone(read);
        if (!cigar.equals(read.getCigar())) {
        	copy.setCigar(cigar);
        	copy.setAttribute(SAMTag.OC.name(), read.getCigarString());
        }
        if (alignmentStart != read.getAlignmentStart()) {
        	copy.setAlignmentStart(alignmentStart);
        	copy.setAttribute(SAMTag.OP.name(), read.getAlignmentStart());
        }
        return copy;
	}
	/**
	 * 0-1 scaled percentage identity of mapped read bases.
	 * @return portion of reference-aligned bases that match reference 
	 */
	public static float getAlignedIdentity(SAMRecord record) {
        Integer nm = record.getIntegerAttribute(SAMTag.NM.name());
		if (nm != null) {
			int refBasesToConsider = record.getReadLength() - SAMRecordUtil.getStartSoftClipLength(record) - SAMRecordUtil.getEndSoftClipLength(record); 
			int refBaseMatches = refBasesToConsider - nm + SequenceUtil.countInsertedBases(record) + SequenceUtil.countDeletedBases(record); 
			return refBaseMatches / (float)refBasesToConsider;
		}
		String md = record.getStringAttribute(SAMTag.MD.name());
		if (StringUtils.isNotEmpty(md)) {
			// Socrates handles this: should we? Which aligners write MD but not NM?
			throw new RuntimeException("Not Yet Implemented: calculation from reads with MD tag but not NM tag");
		}
		throw new IllegalStateException(String.format("Read %s missing NM tag", record.getReadName()));
	}
}
