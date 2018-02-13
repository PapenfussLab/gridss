/**
 * 
 */
package au.edu.wehi.idsv.sam;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.stream.Collectors;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;
import com.google.common.collect.Sets;
import com.google.common.primitives.Booleans;
import com.google.common.primitives.Bytes;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.alignment.AlignerFactory;
import au.edu.wehi.idsv.alignment.Alignment;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.util.IntervalUtil;
import au.edu.wehi.idsv.util.MathUtil;
import au.edu.wehi.idsv.util.MessageThrottler;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamPairUtil;
import htsjdk.samtools.SamPairUtil.PairOrientation;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;

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

	public static int getStartClipLength(SAMRecord r) {
		Cigar cigar = r.getCigar();
		if (cigar == null)
			return 0;
		return getStartClipLength(cigar.getCigarElements());
	}

	public static int getStartClipLength(List<CigarElement> elements) {
		if (elements == null)
			return 0;
		int clipLength = 0;
		for (int i = 0; i < elements.size() && elements.get(i).getOperator().isClipping(); i++) {
			clipLength += elements.get(i).getLength();
		}
		return clipLength;
	}

	public static int getStartSoftClipLength(List<CigarElement> elements) {
		if (elements == null)
			return 0;
		int i = 0;
		while (i < elements.size() && (elements.get(i).getOperator() == CigarOperator.HARD_CLIP
				|| elements.get(i).getOperator() == CigarOperator.SOFT_CLIP)) {
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
		while (i > 0 && (elements.get(i).getOperator() == CigarOperator.HARD_CLIP
				|| elements.get(i).getOperator() == CigarOperator.SOFT_CLIP)) {
			if (elements.get(i).getOperator() == CigarOperator.SOFT_CLIP)
				return elements.get(i).getLength();
			i--;
		}
		return 0;
	}

	public static int getEndClipLength(SAMRecord aln) {
		Cigar cigar = aln.getCigar();
		if (cigar == null)
			return 0;
		return getEndClipLength(cigar.getCigarElements());
	}

	public static int getEndClipLength(List<CigarElement> elements) {
		if (elements == null)
			return 0;
		int clipLength = 0;
		for (int i = elements.size() - 1; i < elements.size() && elements.get(i).getOperator().isClipping(); i--) {
			clipLength += elements.get(i).getLength();
		}
		return clipLength;
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
		seq = Arrays.copyOfRange(seq, record.getReadLength() - getEndSoftClipLength(record), record.getReadLength());
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
		seq = Arrays.copyOfRange(seq, record.getReadLength() - getEndSoftClipLength(record), record.getReadLength());
		return seq;
	}

	public static SAMRecord ensureNmTag(ReferenceSequenceFile ref, SAMRecord record) {
		if (record == null)
			return record;
		if (record.getReadBases() == null)
			return record;
		if (record.getReadBases() == SAMRecord.NULL_SEQUENCE)
			return record;
		if (record.getIntegerAttribute(SAMTag.NM.name()) != null)
			return record;
		if (record.getReadUnmappedFlag())
			return record;
		byte[] refSeq = ref
				.getSubsequenceAt(record.getReferenceName(), record.getAlignmentStart(), record.getAlignmentEnd())
				.getBases();
		final int actualNucleotideDiffs = SequenceUtil.calculateSamNmTag(record, refSeq,
				record.getAlignmentStart() - 1);
		record.setAttribute(SAMTag.NM.name(), actualNucleotideDiffs);
		return record;
	}

	/**
	 * Dovetailing reads either either due to an SV or failure to trim adapters
	 * from a fragment smaller than the read length
	 * 
	 * =======> <=======
	 * 
	 * @param expectedOrientation
	 *            read pair orientation
	 * @return true if the soft clip is due to a fragment size smaller than the
	 *         read length
	 */
	public static boolean isDovetailing(SAMRecord record1, SAMRecord record2, PairOrientation expectedOrientation,
			int margin) {
		if (record1.getReadUnmappedFlag() || record2.getReadUnmappedFlag())
			return false;
		return isDovetailing(record1.getReferenceIndex(), record1.getAlignmentStart(),
				record1.getReadNegativeStrandFlag(), record1.getCigar(), record2.getReferenceIndex(),
				record2.getAlignmentStart(), record2.getReadNegativeStrandFlag(), record2.getCigar(),
				expectedOrientation, margin);
	}

	public static boolean isDovetailing(SAMRecord record, PairOrientation expectedOrientation, int margin) {
		if (record.getReadUnmappedFlag())
			return false;
		if (!record.getReadPairedFlag())
			return false;
		if (record.getMateUnmappedFlag())
			return false;
		Cigar cigar2 = null;
		String mc = record.getStringAttribute(SAMTag.MC.name());
		if (mc != null) {
			cigar2 = TextCigarCodec.decode(mc);
		}
		return isDovetailing(record.getReferenceIndex(), record.getAlignmentStart(), record.getReadNegativeStrandFlag(),
				record.getCigar(), record.getMateReferenceIndex(), record.getMateAlignmentStart(),
				record.getMateNegativeStrandFlag(), cigar2, expectedOrientation, margin);
	}

	private static boolean isDovetailing(int reference1, int start1, boolean isNegativeStrand1, Cigar cigar1,
			int reference2, int start2, boolean isNegativeStrand2, Cigar cigar2, PairOrientation expectedOrientation,
			int margin) {
		if (expectedOrientation != PairOrientation.FR)
			throw new RuntimeException("NYI");
		if (reference1 != reference2)
			return false;
		if (Math.abs(start1 - start2) > margin)
			return false;
		if (isNegativeStrand1 == isNegativeStrand2)
			return false; // FR
		if (cigar2 != null) {
			int end1 = start1 + CigarUtil.referenceLength(cigar1.getCigarElements()) - 1;
			int end2 = start2 + CigarUtil.referenceLength(cigar2.getCigarElements()) - 1;
			if (Math.abs(end1 - end2) > margin)
				return false;
			if (!IntervalUtil.overlapsClosed(start1, end1, start2, end2))
				return false;
		}
		// expect dovetail to look like
		// >>>SSS
		// SSS<<<
		// not like:
		// SSS>>>
		// <<<SSS
		int unexpectedClipLength = isNegativeStrand1 ? getEndSoftClipLength(cigar1.getCigarElements())
				: getStartSoftClipLength(cigar1.getCigarElements());
		if (cigar2 != null) {
			unexpectedClipLength += isNegativeStrand2 ? getEndSoftClipLength(cigar2.getCigarElements())
					: getStartSoftClipLength(cigar2.getCigarElements());
		}
		if (unexpectedClipLength > margin)
			return false;
		return true;
	}

	public static boolean overlap(SAMRecord r1, SAMRecord r2) {
		boolean result = r1 != null && r2 != null && !r1.getReadUnmappedFlag() && !r2.getReadUnmappedFlag()
				&& r1.getReferenceIndex().equals(r2.getReferenceIndex())
				&& ((r1.getAlignmentStart() >= r2.getAlignmentStart() && r1.getAlignmentStart() <= r2.getAlignmentEnd())
						|| (r2.getAlignmentStart() >= r1.getAlignmentStart()
								&& r2.getAlignmentStart() <= r1.getAlignmentEnd()));
		return result;
	}

	/**
	 * Conservative estimate as to whether the reads overlap due to small
	 * fragment size
	 * 
	 * @param local
	 * @param remote
	 * @return number possible breakpoints between the read pair mapped in the
	 *         expected orientation, or Integer.MAX_VALUE if placement is not as
	 *         expected
	 */
	public static boolean estimatedReadsOverlap(SAMRecord record, PairOrientation expectedOrientation,
			int minExpectedAnchorBases) {
		if (expectedOrientation != PairOrientation.FR)
			throw new RuntimeException("NYI");
		if (record.getReadUnmappedFlag() || record.getMateUnmappedFlag()
				|| !record.getReferenceIndex().equals(record.getMateReferenceIndex()) ||
				// FR assumption
				record.getReadNegativeStrandFlag() == record.getMateNegativeStrandFlag()) {
			return false;
		}
		// Assuming FR orientation, adapter sequences have been removed, and no
		// SCs
		if (record.getReadNegativeStrandFlag()) {
			// <----- record
			// -----> mate
			if (record.getMateAlignmentStart() > record.getAlignmentStart())
				return false; // not FR
			// we don't know where the mate actually ends
			// we can't assume = read length as that will
			// remove our pair when it spans a deletion
			// and our mate has an inner soft clip or non-reference CIGAR
			// so we make a conservative guess based on the min bases aligners
			// will return a alignment for
			int offset = Math.min(minExpectedAnchorBases, record.getReadLength()) - 1;
			return record.getMateAlignmentStart() + offset >= record.getAlignmentStart()
					&& record.getMateAlignmentStart() + offset <= record.getAlignmentEnd();

		} else {
			// record ----->
			// <-----
			return record.getMateAlignmentStart() >= record.getAlignmentStart()
					&& record.getMateAlignmentStart() <= record.getAlignmentEnd();
		}
	}

	/**
	 * Estimates the size of sequenced fragment
	 * 
	 * @param record
	 * @return
	 */
	public static int calculateFragmentSize(SAMRecord record1, SAMRecord record2, PairOrientation expectedOrientation) {
		if (expectedOrientation != PairOrientation.FR)
			throw new RuntimeException("NYI");
		if (record1.getReadUnmappedFlag() || record2.getReadUnmappedFlag()
				|| !record1.getReferenceIndex().equals(record2.getReferenceIndex()) ||
				// FR assumption
				record1.getReadNegativeStrandFlag() == record2.getReadNegativeStrandFlag()) {
			return 0;
		}
		// Assuming FR orientation and adapter sequences have been removed
		if (record1.getReadNegativeStrandFlag()) {
			// <--record
			int r1end = record1.getUnclippedEnd();
			int r2start = record2.getUnclippedStart();
			return r1end - r2start + 1;
		} else {
			int r1start = record1.getUnclippedStart();
			int r2end = record2.getUnclippedEnd();
			return r2end - r1start + 1;
		}
	}
	
	/**
	 * Estimates the size of sequenced fragment
	 * 
	 * @param record
	 * @param expectedOrientation
	 * @return
	 */
	public static int estimateFragmentSize(SAMRecord record, PairOrientation expectedOrientation) {
		if (expectedOrientation != PairOrientation.FR)
			throw new RuntimeException("NYI");
		if (record.getReadUnmappedFlag() || record.getMateUnmappedFlag()
				|| !record.getReferenceIndex().equals(record.getMateReferenceIndex()) ||
				// FR assumption
				record.getReadNegativeStrandFlag() == record.getMateNegativeStrandFlag()) {
			return 0;
		}
		// Assuming FR orientation, adapter sequences have been removed
		Cigar mc = SAMUtils.getMateCigar(record);
		if (record.getReadNegativeStrandFlag()) {
			// <--record
			int r1end = record.getUnclippedEnd();
			int r2start = mc == null ?
					// if we don't have a mate cigar we'll just assume that there are is no clipping in the alignment
					record.getMateAlignmentStart() :
					SAMUtils.getMateUnclippedStart(record);
			return r1end - r2start + 1;
		} else {
			int r1start = record.getUnclippedStart();
			int r2end = mc == null ?
					// no MC tag: we have to assume no clipping and the reads are the same length
					record.getMateAlignmentStart() + CigarUtil.readLength(record.getCigar().getCigarElements()) + CigarUtil.countBases(record.getCigar(), CigarOperator.HARD_CLIP) - 1:
					SAMUtils.getMateUnclippedEnd(record);
			return r2end - r1start + 1;
		}
	}
	
	/**
	 * Updates the two reads to indicate the form a read pair
	 * 
	 * @param first
	 *            first in pair
	 * @param second
	 *            second in pair
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
	 * 
	 * @param r
	 *            record to clone
	 * @return clone of SAMRecord
	 */
	public static SAMRecord clone(SAMRecord r) {
		try {
			return (SAMRecord) r.clone();
		} catch (CloneNotSupportedException e) {
			throw new SAMException(e);
		}
	}

	/**
	 * Determine whether the two reads are the same. This is useful for
	 * filtering out read pairs which, due to sequencing chemistry error, both
	 * sequence the same end of the fragment.
	 * 
	 * @param r1
	 *            read to compare
	 * @param r2
	 *            read to compare
	 * @return true if the reads appear to be the same, false otherwise
	 */
	public static boolean areSameRead(SAMRecord r1, SAMRecord r2) {
		// TODO: compare read sequences
		return !r1.getReadUnmappedFlag() && !r2.getReadUnmappedFlag()
				&& r1.getReferenceIndex().equals(r2.getReferenceIndex())
				&& (r1.getUnclippedStart() == r2.getUnclippedStart() || r1.getUnclippedEnd() == r2.getUnclippedEnd())
				&& r1.getReadNegativeStrandFlag() == r2.getReadNegativeStrandFlag();
	}

	/**
	 * Removes soft-clipped bases from the ends of read
	 */
	public static void trimSoftClips(SAMRecord read, int startCount, int endCount) {
		assert (read.getReadUnmappedFlag() || read.getReadLength() == read.getCigar().getReadLength());
		List<CigarElement> cigar = Lists.newArrayList(read.getCigar().getCigarElements());
		if (startCount > 0) {
			CigarElement sc = cigar.get(0);
			int len = sc.getLength();
			assert (sc.getOperator() == CigarOperator.SOFT_CLIP);
			assert (len >= startCount);
			cigar.remove(0);
			len -= startCount;
			if (len > 0) {
				cigar.add(0, new CigarElement(len, CigarOperator.SOFT_CLIP));
			}
		}
		if (endCount > 0) {
			CigarElement sc = cigar.get(cigar.size() - 1);
			int len = sc.getLength();
			assert (sc.getOperator() == CigarOperator.SOFT_CLIP);
			assert (len >= endCount);
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
		assert (read.getReadLength() == read.getCigar().getReadLength());
	}

	public static void trim(SAMRecord read, int startCount, int endCount) {
		int readLength = read.getReadLength();
		read.setReadBases(Arrays.copyOfRange(read.getReadBases(), startCount, readLength - endCount));
		read.setBaseQualities(Arrays.copyOfRange(read.getBaseQualities(), startCount, readLength - endCount));
		if (read.getCigar() != null && read.getCigar().getCigarElements().size() > 0) {
			read.setAlignmentStart(read.getAlignmentStart() + CigarUtil.offsetOf(read.getCigar(), startCount));
			read.setCigar(CigarUtil.trimReadBases(read.getCigar(), startCount, endCount));
			assert (read.getReadLength() == read.getCigar().getReadLength());
		}
	}

	public static boolean isReferenceAlignment(SAMRecord r) {
		if (r.getCigar() == null)
			return true;
		return isReferenceAlignment(r.getCigar().getCigarElements());
	}

	public static boolean isReferenceAlignment(List<CigarElement> cigar) {
		if (cigar == null || cigar.size() == 0)
			return true;
		for (CigarElement e : cigar) {
			if (e.getLength() > 0) {
				switch (e.getOperator()) {
				case M:
				case EQ:
				case X:
				case P:
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
	 * 
	 * @param read
	 *            read
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
	 * 
	 * @param read
	 *            read
	 * @return number of bits of entropy
	 */
	public static double entropy(SAMRecord read) {
		double bits = au.edu.wehi.idsv.util.SequenceUtil.shannonEntropy(read.getReadBases(), 0, read.getReadLength());
		return bits;
	}

	/**
	 * Performs local realignment of the given SAMRecord in a window around the
	 * alignment location
	 * 
	 * @param reference
	 *            reference genome
	 * @param read
	 *            read
	 * @param windowSize
	 *            number of bases to extend window around alignment
	 * @param extendWindowForClippedBases
	 *            extend window for soft clipped bases
	 * @return copy of realigned read if alignment changed, read if realignment
	 *         did not change the alignment
	 */
	public static SAMRecord realign(ReferenceLookup reference, SAMRecord read, int windowSize,
			boolean extendWindowForClippedBases) {
		if (read.getReadUnmappedFlag())
			return read;
		SAMSequenceRecord refSeq = reference.getSequenceDictionary().getSequence(read.getReferenceIndex());
		// find reference bounds of read. Negative deletions mean we can't just
		// use
		// getAlignmentStart() and getAlignmentEnd()
		int pos = read.getAlignmentStart();
		int start = pos;
		int end = pos;
		for (CigarElement ce : read.getCigar().getCigarElements()) {
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
		// defensive checks so we don't crash the JVM if an unexpected character
		// is encountered
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
		Alignment alignment = null;
		if (ass != null && ass.length > 0 && ref != null && ref.length > 0) {
			try {
				alignment = AlignerFactory.create().align_smith_waterman(ass, ref);
			} catch (Exception e) {
				// swallow and log alignment error
				if (!MessageThrottler.Current.shouldSupress(log, "local realignment failures")) {
					log.error(e, String.format("Error aligning %s to %s:%d-%d", read.getReadName(),
							refSeq.getSequenceName(), start, end));
				}
			}
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
	 * 
	 * @return portion of reference-aligned bases that match reference
	 */
	public static float getAlignedIdentity(SAMRecord record) {
		Integer nm = record.getIntegerAttribute(SAMTag.NM.name());
		if (nm != null) {
			int refBasesToConsider = record.getReadLength() - SAMRecordUtil.getStartSoftClipLength(record)
					- SAMRecordUtil.getEndSoftClipLength(record);
			int refBaseMatches = refBasesToConsider - nm + SequenceUtil.countInsertedBases(record)
					+ SequenceUtil.countDeletedBases(record);
			return refBaseMatches / (float) refBasesToConsider;
		}
		String md = record.getStringAttribute(SAMTag.MD.name());
		if (StringUtils.isNotEmpty(md)) {
			// Socrates handles this: should we? Which aligners write MD but not
			// NM?
			throw new RuntimeException("Not Yet Implemented: calculation from reads with MD tag but not NM tag");
		}
		throw new IllegalStateException(String.format("Read %s missing NM tag", record.getReadName()));
	}

	/**
	 * Computed tags requiring reads from the same template to be processed
	 * together
	 */
	public static final Set<String> TEMPLATE_TAGS = ImmutableSet.of(SAMTag.CC.name(), SAMTag.CP.name(), SAMTag.FI.name(), SAMTag.HI.name(),
			SAMTag.IH.name(), SAMTag.Q2.name(), SAMTag.R2.name(), SAMTag.MC.name(), SAMTag.MQ.name(), SAMTag.SA.name(), SAMTag.TC.name());

	/**
	 * Populates the missing template tags for the records originating from the
	 * same template
	 * 
	 * @param records
	 *            Records from the same template.
	 * @param tags
	 *            Tags to populate
	 * @param restoreHardClips
	 *            convert hard clips into soft clips if the sequence is
	 *            available from a chimeric alignment of the same segmentt
	 * @param fixMates
	 * @param recalculateSupplementary 
	 */
	public static void calculateTemplateTags(List<SAMRecord> records, Set<String> tags,
			boolean restoreHardClips, boolean fixMates, boolean recalculateSupplementary) {
		List<List<SAMRecord>> segments = templateBySegment(records);
		// FI
		if (tags.contains(SAMTag.FI.name())) {
			for (int i = 0; i < segments.size(); i++) {
				for (SAMRecord r : segments.get(i)) {
					r.setAttribute(SAMTag.FI.name(), i);
				}
			}
		}
		// TC
		if (tags.contains(SAMTag.TC.name())) {
			for (SAMRecord r : records) {
				r.setAttribute(SAMTag.TC.name(), segments.size());
			}
		}
		if (restoreHardClips) {
			for (int i = 0; i < segments.size(); i++) {
				softenHardClips(segments.get(i));
			}
		}
		// SA
		if (tags.contains(SAMTag.SA.name())) {
			for (int i = 0; i < segments.size(); i++) {
				calculateSATags(segments.get(i));
			}
		}
		if (recalculateSupplementary) {
			for (int i = 0; i < segments.size(); i++) {
				recalculateSupplementaryFromSA(segments.get(i));
			}
		}
		// R2
		if (tags.contains(SAMTag.R2.name())) {
			for (int i = 0; i < segments.size(); i++) {
				byte[] br2 = getConsensusSequence(segments.get((i + 1) % segments.size()));
				String r2 = br2 != null && br2 != SAMRecord.NULL_SEQUENCE &&  segments.size() > 1 ? StringUtil.bytesToString(br2) : null;
				for (SAMRecord r : segments.get(i)) {
					r.setAttribute(SAMTag.R2.name(), r2);
				}
			}
		}
		// Q2
		if (tags.contains(SAMTag.Q2.name())) {
			for (int i = 0; i < segments.size(); i++) {
				byte[] bq2 = getConsensusBaseQualities(segments.get((i + 1) % segments.size()));
				String q2 = bq2 != null && bq2 != SAMRecord.NULL_QUALS && segments.size() > 1 ? SAMUtils.phredToFastq(bq2) : null;
				for (SAMRecord r : segments.get(i)) {
					r.setAttribute(SAMTag.Q2.name(), q2);
				}
			}
		}
		if (Sets.intersection(tags,
				ImmutableSet.of(SAMTag.CC.name(), SAMTag.CP.name(), SAMTag.HI.name(), SAMTag.IH.name())).size() > 0) {
			for (int i = 0; i < segments.size(); i++) {
				calculateMultimappingTags(tags, segments.get(i));
			}
		}
		if (fixMates || tags.contains(SAMTag.MC.name()) || tags.contains(SAMTag.MQ.name())) {
			if (segments.size() >= 2) {
				// hack to prevent SAMRecord getters from throwing an exception
				records.stream().forEach(r -> r.setReadPairedFlag(true));
				segments.get(0).stream().forEach(r -> r.setFirstOfPairFlag(true));
				segments.get(1).stream().forEach(r -> r.setSecondOfPairFlag(true));
				fixMates(segments, tags.contains(SAMTag.MC.name()));
			} else {
				for (SAMRecord r : records) {
					clearMateInformation(r, true);
				}
			}
		}
		if (tags.contains(SamTags.MULTIMAPPING_FRAGMENT)) {
			int mappingLocations = 0;
			for (SAMRecord r : records) {
				if (!r.getSupplementaryAlignmentFlag()) {
					mappingLocations++;
				}
			}
			for (SAMRecord r : records) {
				r.setAttribute(SamTags.MULTIMAPPING_FRAGMENT, mappingLocations <= segments.size() ? null : mappingLocations - segments.size());
			}
		}
	}
	/**
	 * Orders the records such that the primary record for a split read alignment is first
	 */
	private static Ordering<SAMRecord> ByBestPrimarySplitCandidate = new Ordering<SAMRecord>() {
		public int compare(SAMRecord arg1, SAMRecord arg2) {
			return ComparisonChain.start()
				// already flagged as supp is bad  
				.compareFalseFirst(arg1.getSupplementaryAlignmentFlag(), arg2.getSupplementaryAlignmentFlag())
				// flagged as secondary is bad due to legacy treatment of secondary alignments as supplementary (eg bwa mem -M) 
				.compareFalseFirst(arg1.isSecondaryAlignment(), arg2.isSecondaryAlignment()) 
				// the record with the shorter soft clip is a better candidate
				.compare(SAMRecordUtil.getStartClipLength(arg1) + SAMRecordUtil.getEndClipLength(arg1), SAMRecordUtil.getStartClipLength(arg2) + SAMRecordUtil.getEndClipLength(arg2))
				// Other options are:
				// - record is flagged as the mate of a read pair (strong support for that record to be the primary)
				// - best MAPQ
				.result();
		}
	};
	private static void recalculateSupplementaryFromSA(List<SAMRecord> segments) {
		HashMap<List<ChimericAlignment>, List<SAMRecord>> saLookup = new HashMap<>();
		for (SAMRecord r : segments) {
			List<ChimericAlignment> splitca = ChimericAlignment.getChimericAlignments(r);
			if (splitca.isEmpty() || r.getReadUnmappedFlag()) {
				r.setSupplementaryAlignmentFlag(false);
			} else {
				splitca.add(new ChimericAlignment(r));
				splitca.sort(ChimericAlignment.ByReadOffset);
				List<SAMRecord> saGroup = saLookup.get(splitca);
				if (saGroup == null) {
					saGroup = new ArrayList<>();
					saLookup.put(splitca, saGroup);
				}
				saGroup.add(r);
			}
		}
		for (List<SAMRecord> split : saLookup.values()) {
			split.sort(ByBestPrimarySplitCandidate);
			for (int i = 0; i < split.size(); i++) {
				split.get(i).setSupplementaryAlignmentFlag(i != 0);
			}
		}
	}

	private static SAMRecord getPrimarySplitAlignmentFor(SAMRecord rec, List<SAMRecord> options) {
		List<ChimericAlignment> splits = ChimericAlignment.getChimericAlignments(rec);
		// supposed to be the first split alignment
		final ChimericAlignment primary = splits.size() == 0 ? null : splits.get(0);
		Optional<SAMRecord> best = options.stream()
				.filter(r -> !r.getSupplementaryAlignmentFlag())
				.filter(r -> new ChimericAlignment(rec).equals(primary))
				.findFirst();
		// try any split alignment
		if (!best.isPresent()) {
			best = options.stream()
					.filter(r -> !r.getSupplementaryAlignmentFlag())
					.filter(r -> splits.contains(new ChimericAlignment(rec)))
					.findFirst();
		}
		// just grab the primary
		if (!best.isPresent()) {
			best = options.stream()
					.filter(r -> !r.getSupplementaryAlignmentFlag())
					.filter(r -> r.isSecondaryAlignment() == rec.isSecondaryAlignment())
					.findFirst();
		}
		// grab anything that's not a supplementary
		if (!best.isPresent()) {
			best = options.stream()
					.filter(r -> !r.getSupplementaryAlignmentFlag())
					.findFirst();
		}
		// Just use ourself if there's nothing else out there but other
		// supplementary alignments
		return best.orElse(rec);
	}

	private static void fixMates(List<List<SAMRecord>> segments, boolean mateCigar) {
		for (int i = 0; i < segments.size(); i++) {
			List<SAMRecord> currentSegment = segments.get(i);
			List<SAMRecord> nextSegment = segments.get((i + 1) % segments.size());
			if (segments.size() == 1) {
				// don't set ourselves as a mate
				nextSegment = Collections.emptyList();
			}
			// resort so we match up the primary records last
			currentSegment.sort(Comparator.comparing(SAMRecord::isSecondaryAlignment).reversed());
			for (SAMRecord r : currentSegment) {
				if (r.getSupplementaryAlignmentFlag()) {
					SAMRecord primary = getPrimarySplitAlignmentFor(r, currentSegment);
					SAMRecord mate = bestMateFor(primary, nextSegment);
					if (mate != null) {
						SamPairUtil.setMateInformationOnSupplementalAlignment(r, mate, mateCigar);
					} else {
						clearMateInformation(r, true);
					}
				} else {
					SAMRecord mate = bestMateFor(r, nextSegment);
					if (mate != null) {
						SamPairUtil.setMateInfo(r, mate, mateCigar);
					} else {
						clearMateInformation(r, true);
					}
				}
			}
		}
	}

	/**
	 * Ensures that the read is no longer considered paired
	 * 
	 * @param r
	 *            read to clear mate information from
	 * @param fullClear
	 *            strips all mate information from the read
	 */
	private static void clearMateInformation(SAMRecord r, boolean fullClear) {
		r.setReadPairedFlag(false);
		if (fullClear) {
			r.setReadPairedFlag(true);
			if (r.getFirstOfPairFlag()) {
				r.setAttribute(SAMTag.FI.name(), 0);
			} else if (r.getSecondOfPairFlag()) {
				r.setAttribute(SAMTag.FI.name(), 1);
			}
			r.setReadPairedFlag(false);
			r.setFirstOfPairFlag(false);
			r.setSecondOfPairFlag(false);
			r.setProperPairFlag(false);
			r.setMateUnmappedFlag(false);
			r.setMateNegativeStrandFlag(false);
			r.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
			r.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
			r.setAttribute(SAMTag.MC.name(), null);
			r.setAttribute(SAMTag.MQ.name(), null);
			r.setAttribute(SAMTag.R2.name(), null);
			r.setAttribute(SAMTag.Q2.name(), null);
		}
	}

	private static SAMRecord bestMateFor(SAMRecord r, List<SAMRecord> mates) {
		// see if we already point to one of the mates (ie: don't touch existing records)
		Optional<SAMRecord> best = mates.stream().filter(m -> m.getReadUnmappedFlag() == r.getMateUnmappedFlag())
				.filter(m -> m.getReferenceIndex() == r.getMateReferenceIndex())
				.filter(m -> m.getAlignmentStart() == r.getMateAlignmentStart())
				.filter(m -> m.getReadNegativeStrandFlag() == r.getMateNegativeStrandFlag()).findFirst();
		// see if one of the mates already points to us
		if (best.isPresent()) return best.get();
		best = mates.stream().filter(m -> r.getReadUnmappedFlag() == m.getMateUnmappedFlag())
				.filter(m -> r.getReferenceIndex() == m.getMateReferenceIndex())
				.filter(m -> r.getAlignmentStart() == m.getMateAlignmentStart())
				.filter(m -> r.getReadNegativeStrandFlag() == m.getMateNegativeStrandFlag()).findFirst();
		if (best.isPresent()) return best.get();
		// match primary with primary and secondary with secondary
		best = mates.stream().filter(m -> m.getSupplementaryAlignmentFlag() == r.getSupplementaryAlignmentFlag())
				.filter(m -> m.isSecondaryAlignment() == r.isSecondaryAlignment())
				// prefer mapped reads
				.sorted(Comparator.comparing(SAMRecord::getReadUnmappedFlag)).findFirst();
		if (best.isPresent()) return best.get();
		// just get anything that's not a supplementary alignment
		best = mates.stream().filter(m -> m.getSupplementaryAlignmentFlag() == r.getSupplementaryAlignmentFlag())
				.findFirst();
		if (best.isPresent()) return best.get();
		// we give up - anything will do at this point
		best = mates.stream().findFirst();
		return best.orElse(null);
	}

	private static final List<List<SAMRecord>> templateBySegment(List<SAMRecord> records) {
		List<List<SAMRecord>> segments = new ArrayList<List<SAMRecord>>();
		for (SAMRecord r : records) {
			int index = SAMRecordUtil.getSegmentIndex(r);
			while (index >= segments.size()) {
				segments.add(new ArrayList<SAMRecord>());
			}
			segments.get(index).add(r);
		}
		return segments;
	}

	public static void calculateMultimappingTags(Set<String> tags, List<SAMRecord> list) {
		// TODO: how does SA split read alignment interact with multimapping
		// reads?
		list.sort(new SAMRecordCoordinateComparator());
		for (int i = 0; i < list.size(); i++) {
			SAMRecord r = list.get(i);
			if (tags.contains(SAMTag.IH.name())) {
				// TODO: how does this differ from NH? Does NH count unreported
				// alignments?
				r.setAttribute(SAMTag.IH.name(), list.size());
			}
			if (tags.contains(SAMTag.HI.name())) {
				r.setAttribute(SAMTag.HI.name(), i);
			}
			if (tags.contains(SAMTag.CC.name())) {
				r.setAttribute(SAMTag.CC.name(), list.get((i + 1) % list.size()).getReferenceName());
			}
			if (tags.contains(SAMTag.CP.name())) {
				r.setAttribute(SAMTag.CP.name(), list.get((i + 1) % list.size()).getAlignmentStart());
			}
		}
	}

	private static void calculateSATags(List<SAMRecord> list) {
		// TODO use CC, CP, HI, IH tags if they are present
		// TODO break ties in an other other than just overwriting the previous
		// alignment that
		// finishes at a given position
		SortedMap<Integer, SAMRecord> start = new TreeMap<Integer, SAMRecord>();
		SortedMap<Integer, SAMRecord> end = new TreeMap<Integer, SAMRecord>();
		for (SAMRecord r : list) {
			start.put(getFirstAlignedBaseReadOffset(r), r);
			end.put(getLastAlignedBaseReadOffset(r), r);
		}
		for (SAMRecord r : list) {
			if (r.getAttribute(SAMTag.SA.name()) != null)
				continue;
			ArrayDeque<SAMRecord> alignment = new ArrayDeque<SAMRecord>();
			SAMRecord prev = r;
			while (prev != null) {
				SortedMap<Integer, SAMRecord> before = end.headMap(getFirstAlignedBaseReadOffset(prev));
				if (before.isEmpty()) {
					prev = null;
				} else {
					prev = before.get(before.lastKey());
					alignment.addFirst(prev);
				}
			}
			SAMRecord next = r;
			while (next != null) {
				SortedMap<Integer, SAMRecord> after = start.tailMap(getLastAlignedBaseReadOffset(next) + 1);
				if (after.isEmpty()) {
					next = null;
				} else {
					next = after.get(after.firstKey());
					alignment.addLast(next);
				}
			}
			if (alignment.size() == 0) {
				r.setAttribute(SAMTag.SA.name(), null);
			} else {
				List<SAMRecord> aln = new ArrayList<>(alignment);
				Collections.sort(aln, new Ordering<SAMRecord>() {
					@Override
					public int compare(SAMRecord arg0, SAMRecord arg1) {
						// http://samtools.github.io/hts-specs/SAMtags.pdf
						// Conventionally, at a supplementary line, the first
						// element points to the primary line.
						return Booleans.compare(arg0.getSupplementaryAlignmentFlag(),
								arg1.getSupplementaryAlignmentFlag());
					}
				});
				StringBuilder sb = new StringBuilder();
				for (SAMRecord chim : aln) {
					sb.append(new ChimericAlignment(chim));
					sb.append(";");
				}
				sb.deleteCharAt(sb.length() - 1);
				r.setAttribute(SAMTag.SA.name(), sb.toString());
			}
		}
		Set<ChimericAlignment> referencedReads = list.stream()
				.flatMap(r -> ChimericAlignment.getChimericAlignments(r.getStringAttribute(SAMTag.SA.name())).stream())
				.collect(Collectors.toSet());
		// validate SA tags
		for (SAMRecord r : list) {
			referencedReads.remove(new ChimericAlignment(r));
		}
		if (!referencedReads.isEmpty()) {
			if (!MessageThrottler.Current.shouldSupress(log, "Missing records referred to by SA tags")) {
				log.warn(String.format("SA tag of read %s refers to missing alignments %s", list.get(0).getReadName(), referencedReads));
			}
		}
	}

	/**
	 * Converts any hard clips into soft clips using alternate alignments
	 * 
	 * @param records
	 */
	public static final void softenHardClips(List<SAMRecord> records) {
		for (SAMRecord r : records) {
			hardClipToN(r);
		}
		byte[] seq = getConsensusSequence(records);
		byte[] qual = getConsensusBaseQualities(records);
		for (SAMRecord r : records) {
			Cigar c = r.getCigar();
			if (r.getReadUnmappedFlag() || c.getReadLength() == seq.length) {
				if (seq != null && seq != SAMRecord.NULL_SEQUENCE) {
					byte[] newseq = Arrays.copyOf(seq, seq.length);
					if (r.getReadNegativeStrandFlag()) {
						SequenceUtil.reverseComplement(newseq);
					}
					r.setReadBases(newseq);
				}
				if (qual != null && qual != SAMRecord.NULL_QUALS) {
					byte[] newqual = Arrays.copyOf(qual, qual.length);
					if (r.getReadNegativeStrandFlag()) {
						ArrayUtils.reverse(newqual);
					}
					r.setBaseQualities(newqual);
				}
			} else {
				if (!MessageThrottler.Current.shouldSupress(log, "softening hard clips")) {
					log.warn(String.format("Input sanity check failure: different alignment records imply different read lengths %s. "
							+ "This can be cause by GATK indel realignment incorrectly removing hard clipped bases when realigning. Do not use GATK"
							+ " indel realigned BAM files with GRIDSS.", r.getReadName()));
				}
			}
		}
	}

	/**
	 * Converts any hard clips into soft clipped N bases
	 * 
	 * @param record
	 */
	public static final void hardClipToN(SAMRecord r) {
		if (r.getReadUnmappedFlag() || r.getCigar() == null)
			return;
		if (!Iterables.any(r.getCigar().getCigarElements(), ce -> ce.getOperator() == CigarOperator.HARD_CLIP))
			return;
		List<CigarElement> list = new ArrayList<>(r.getCigar().getCigarElements());
		int startlength = r.getCigar().getFirstCigarElement().getOperator() == CigarOperator.HARD_CLIP
				? r.getCigar().getFirstCigarElement().getLength() : 0;
		int endlength = r.getCigar().getLastCigarElement().getOperator() == CigarOperator.HARD_CLIP
				? r.getCigar().getLastCigarElement().getLength() : 0;
		r.setCigar(new Cigar(CigarUtil.clean(list.stream()
				.map(ce -> ce.getOperator() == CigarOperator.HARD_CLIP
						? new CigarElement(ce.getLength(), CigarOperator.SOFT_CLIP) : ce)
				.collect(Collectors.toList()))));
		if (r.getReadBases() != null && r.getReadBases() != SAMRecord.NULL_SEQUENCE) {
			byte[] start = new byte[startlength];
			byte[] end = new byte[endlength];
			Arrays.fill(start, (byte) 'N');
			Arrays.fill(end, (byte) 'N');
			r.setReadBases(Bytes.concat(start, r.getReadBases(), end));
		}
		if (r.getBaseQualities() != null && r.getBaseQualities() != SAMRecord.NULL_QUALS) {
			byte[] start = new byte[startlength];
			byte[] end = new byte[endlength];
			Arrays.fill(start, (byte) 0);
			Arrays.fill(end, (byte) 0);
			r.setBaseQualities(Bytes.concat(start, r.getBaseQualities(), end));
		}
	}

	/**
	 * Gets the full read sequence
	 * 
	 * @param records
	 * @return
	 */
	private static final byte[] getConsensusSequence(List<SAMRecord> records) {
		if (records == null || records.size() == 0)
			return null;
		byte[] cons = SAMRecord.NULL_SEQUENCE;
		for (SAMRecord r : records) {
			boolean negativeStrand = r.getReadNegativeStrandFlag();
			byte[] seq = r.getReadBases();
			if (seq == null || seq.length == 0) continue;
			if (seq.length > cons.length) {
				cons = Arrays.copyOf(seq, seq.length);
				if (negativeStrand) {
					SequenceUtil.reverseComplement(cons);
				}
			} else {
				for (int i = 0; i < seq.length; i++) {
					if (cons[i] == 'N') {
						int seqOffset = negativeStrand ? seq.length - 1 - i : i;
						cons[i] = negativeStrand ? SequenceUtil.complement(seq[seqOffset]) : seq[seqOffset];
					}
				}
			}
		}
		return cons;
	}
	private static final byte[] getConsensusBaseQualities(List<SAMRecord> records) {
		if (records == null || records.size() == 0)
			return null;
		byte[] cons = SAMRecord.NULL_QUALS;
		for (SAMRecord r : records) {
			boolean negativeStrand = r.getReadNegativeStrandFlag();
			byte[] qual = r.getBaseQualities();
			if (qual == null || qual.length == 0) continue;
			if (qual.length > cons.length) {
				cons = Arrays.copyOf(qual, qual.length);
				if (negativeStrand) {
					ArrayUtils.reverse(cons);
				}
			} else {
				for (int i = 0; i < qual.length; i++) {
					int offset = negativeStrand ? qual.length - 1 - i : i;
					cons[i] = (byte) Math.max(cons[i], qual[offset]);
				}
			}
		}
		return cons;
	}
	/**
	 * The index of segment in the template
	 * 
	 * @param record
	 *            SAMRecord
	 * @return The index of segment in the template
	 */
	public static final int getSegmentIndex(SAMRecord record) {
		Integer hiTag = record.getIntegerAttribute(SAMTag.FI.name());
		if (hiTag != null)
			return hiTag;
		if (record.getReadPairedFlag()) {
			if (record.getFirstOfPairFlag())
				return 0;
			if (record.getSecondOfPairFlag())
				return 1;
		}
		return 0; // default to the first segment as per sam specs
	}

	/**
	 * Gets the read offset of the first aligned base.
	 * 
	 * @param r
	 * @return zero based offset from start of read based on fastq sequencing
	 *         base order
	 */
	public static final int getFirstAlignedBaseReadOffset(SAMRecord r) {
		Cigar c = r.getCigar();
		if (c == null || c.getCigarElements().size() == 0)
			return -1;
		if (r.getReadNegativeStrandFlag()) {
			return getEndClipLength(c.getCigarElements());
		} else {
			return getStartClipLength(c.getCigarElements());
		}
	}

	/**
	 * Gets the read offset of the final aligned base.
	 * 
	 * @param r
	 * @return zero based offset of the final aligned read base based on fastq
	 *         sequencing base order
	 */
	public static final int getLastAlignedBaseReadOffset(SAMRecord r) {
		Cigar c = r.getCigar();
		if (c == null || c.getCigarElements().size() == 0)
			return -1;
		int length = c.getReadLength();
		if (c.getFirstCigarElement().getOperator() == CigarOperator.HARD_CLIP) {
			length += c.getFirstCigarElement().getLength();
		}
		if (c.getLastCigarElement().getOperator() == CigarOperator.HARD_CLIP) {
			length += c.getLastCigarElement().getLength();
		}

		if (r.getReadNegativeStrandFlag()) {
			return length - getStartClipLength(c.getCigarElements()) - 1;
		} else {
			return length - getEndClipLength(c.getCigarElements()) - 1;

		}
	}

	/**
	 * Converts a record with mapq below the given mapq threshold to unmapped
	 * reads
	 * 
	 * @param record
	 *            SAMRecord to convert
	 * @param minMapq
	 *            minimum MAPQ to avoid unmapping
	 */
	public static SAMRecord lowMapqToUnmapped(SAMRecord record, double minMapq) {
		boolean primaryUnmapped = lowMapqToUnmapped_SA(record, minMapq);
		if (record.getMappingQuality() < minMapq || primaryUnmapped) {
			record.setReadUnmappedFlag(true);
		}
		if (record.getReadPairedFlag()) {
			Integer mateMapq = record.getIntegerAttribute(SAMTag.MQ.name());
			if (mateMapq != null && mateMapq < minMapq) {
				record.setMateUnmappedFlag(true);
			}
		}
		return record;
	}

	private static boolean lowMapqToUnmapped_SA(SAMRecord record, double minMapq) {
		if (record.getAttribute(SAMTag.SA.name()) == null) {
			return false;
		}
		
		List<ChimericAlignment> alignments = ChimericAlignment.getChimericAlignments(record);
		// By convention, primary is the first SA record
		boolean primaryUnmapped = record.isSecondaryOrSupplementary() &&
				(alignments.size() == 0 || alignments.get(0).mapq < minMapq);
		record.setAttribute(SAMTag.SA.name(), alignments.stream()
				.filter(ca -> ca.mapq >= minMapq)
				.map(ca -> ca.toString())
				.collect(Collectors.joining(";")));
		return primaryUnmapped;
	}

	/**
	 * Converts soft clipped bases that match the reference genome to alignment
	 * matches. In general, such bases should be included in the alignment but,
	 * due to error correction of sequencing errors, this does not necessarily
	 * hold true for assembly contigs.
	 * 
	 * @param r
	 *            record
	 */
	public static void unclipExactReferenceMatches(ReferenceLookup ref, SAMRecord read) {
		if (read.getReadUnmappedFlag())
			return;
		byte[] seq = read.getReadBases();
		int refIndex = read.getReferenceIndex();
		SAMSequenceRecord refseq = ref.getSequenceDictionary().getSequence(refIndex);
		int startunclip = 0;
		int startclip = getStartSoftClipLength(read);
		for (int i = startclip - 1; i >= 0; i--) {
			int pos = read.getAlignmentStart() - startclip + i;
			if (pos >= 1 && pos <= refseq.getSequenceLength()) {
				byte refbase = ref.getBase(refIndex, pos);
				byte readbase = seq[i];
				if (SequenceUtil.basesEqual(refbase, readbase) && !SequenceUtil.basesEqual(SequenceUtil.N, readbase)) {
					startunclip++;
				} else {
					break;
				}
			}
		}
		int endunclip = 0;
		int endcliplength = getEndSoftClipLength(read);
		for (int i = 0; i < endcliplength; i++) {
			int pos = read.getAlignmentEnd() + i + 1;
			if (pos >= 1 && pos <= refseq.getSequenceLength()) {
				byte refbase = ref.getBase(refIndex, read.getAlignmentEnd() + i + 1);
				byte readbase = seq[read.getReadLength() - endcliplength + i];
				if (SequenceUtil.basesEqual(refbase, readbase)) {
					endunclip++;
				} else {
					break;
				}
			}
		}
		read.setCigar(softClipToAligned(read.getCigar(), startunclip, endunclip));
		read.setAlignmentStart(read.getAlignmentStart() - startunclip);
	}

	private static Cigar softClipToAligned(Cigar cigar, int startunclip, int endunclip) {
		if (startunclip == 0 && endunclip == 0)
			return cigar;
		List<CigarElement> list = cigar.getCigarElements();
		list = CigarUtil.clean(list, true);
		if (startunclip > 0) {
			int i = 0;
			while (i < list.size() && list.get(i).getOperator() != CigarOperator.SOFT_CLIP)
				i++;
			if (i == list.size() || list.get(i).getLength() < startunclip) {
				throw new IllegalArgumentException(
						String.format("Unable to unclip %d clipped bases from start of %s", startunclip, cigar));
			}
			list.set(i, new CigarElement(list.get(i).getLength() - startunclip, CigarOperator.SOFT_CLIP));
			list.add(i + 1, new CigarElement(startunclip, CigarOperator.MATCH_OR_MISMATCH));
			CigarUtil.clean(list, false);
		}
		if (endunclip > 0) {
			int i = list.size() - 1;
			while (i > 0 && list.get(i).getOperator() != CigarOperator.SOFT_CLIP)
				i--;
			if (i < 0 || list.get(i).getLength() < endunclip) {
				throw new IllegalArgumentException(
						String.format("Unable to unclip %d clipped bases from end of %s", endunclip, cigar));
			}
			list.set(i, new CigarElement(list.get(i).getLength() - endunclip, CigarOperator.SOFT_CLIP));
			list.add(i - 1, new CigarElement(endunclip, CigarOperator.MATCH_OR_MISMATCH));
			CigarUtil.clean(list, false);
		}
		return new Cigar(list);
	}
	public static double getEffectiveMapq(SAMRecord r, double fallbackMapq) {
		if (r.getReadUnmappedFlag()) return 0;
		Integer hits = reportedAlignments(r);
		if (hits != null) {
			if (hits > 1) {
				// consider all mappings equally likely since we don't have any more information than that
				return MathUtil.prToPhred((hits - 1) / (double)hits);
			} else {
				// TODO: what do to here?
				// bowtie2 doesn't conform to the specs in multi-mapping mode
				// and reports a nonsense mapq
				return fallbackMapq;
			}
		}
		int mapq = r.getMappingQuality();
		if (mapq == SAMRecord.UNKNOWN_MAPPING_QUALITY) {
			return fallbackMapq;
		}
		return mapq;
	}
	private static Integer reportedAlignments(SAMRecord r) {
		if (r.getIntegerAttribute(SAMTag.NH.name()) != null) {
			// NH i Number of reported alignments that contains the query in the current record
			return r.getIntegerAttribute(SAMTag.NH.name());
		}
		if (r.getIntegerAttribute(SAMTag.IH.name()) != null) {
			// IH i Number of stored alignments in SAM that contains the query in the current record
			return r.getIntegerAttribute(SAMTag.IH.name());
		}
		return null;
	}
}
