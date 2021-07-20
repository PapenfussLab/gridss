/**
 * 
 */
package au.edu.wehi.idsv.sam;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.alignment.AlignerFactory;
import au.edu.wehi.idsv.alignment.Alignment;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.util.GroupingIterator;
import au.edu.wehi.idsv.util.IntervalUtil;
import au.edu.wehi.idsv.util.MathUtil;
import au.edu.wehi.idsv.util.MessageThrottler;
import com.google.common.collect.*;
import com.google.common.collect.ImmutableRangeSet.Builder;
import com.google.common.primitives.Booleans;
import com.google.common.primitives.Bytes;
import com.google.common.primitives.Ints;
import htsjdk.samtools.*;
import htsjdk.samtools.SamPairUtil.PairOrientation;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.util.Pair;

import java.util.*;
import java.util.stream.Collectors;

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
		if (elements == null) {
			return 0;
		}
		int clipLength = 0;
		for (int i = 0; i < elements.size() && elements.get(i).getOperator().isClipping(); i++) {
			clipLength += elements.get(i).getLength();
		}
		return clipLength;
	}

	public static int getStartSoftClipLength(List<CigarElement> elements) {
		if (elements == null) {
			return 0;
		}
		int i = 0;
		while (i < elements.size() && elements.get(i).getOperator().isClipping()) {
			if (elements.get(i).getOperator() == CigarOperator.SOFT_CLIP) {
				return elements.get(i).getLength();
			}
			i++;
		}
		return 0;
	}

	public static int getEndSoftClipLength(SAMRecord aln) {
		Cigar cigar = aln.getCigar();
		if (cigar == null) {
			return 0;
		}
		return getEndSoftClipLength(cigar.getCigarElements());
	}

	public static int getEndSoftClipLength(List<CigarElement> elements) {
		if (elements == null) return 0;
		int i = elements.size() - 1;
		while (i >= 0 && elements.get(i).getOperator().isClipping()) {
			if (elements.get(i).getOperator() == CigarOperator.SOFT_CLIP) {
				return elements.get(i).getLength();
			}
			i--;
		}
		return 0;
	}

	public static int getEndClipLength(SAMRecord aln) {
		Cigar cigar = aln.getCigar();
		if (cigar == null) {
			return 0;
		}
		return getEndClipLength(cigar.getCigarElements());
	}

	public static int getEndClipLength(List<CigarElement> elements) {
		if (elements == null) {
			return 0;
		}
		int clipLength = 0;
		for (int i = elements.size() - 1; i >= 0 && elements.get(i).getOperator().isClipping(); i--) {
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
		if (seq == null) {
			return null;
		}
		if (seq == SAMRecord.NULL_SEQUENCE) {
			return null;
		}
		seq = Arrays.copyOfRange(seq, 0, getStartSoftClipLength(record));
		return seq;
	}

	public static byte[] getEndSoftClipBases(SAMRecord record) {
		byte[] seq = record.getReadBases();
		if (seq == null) {
			return null;
		}
		if (seq == SAMRecord.NULL_SEQUENCE) {
			return null;
		}
		seq = Arrays.copyOfRange(seq, record.getReadLength() - getEndSoftClipLength(record), record.getReadLength());
		return seq;
	}

	public static byte[] getStartSoftClipBaseQualities(SAMRecord record) {
		byte[] seq = record.getBaseQualities();
		if (seq == null) {
			return null;
		}
		if (seq == SAMRecord.NULL_QUALS) {
			return null;
		}
		seq = Arrays.copyOfRange(seq, 0, getStartSoftClipLength(record));
		return seq;
	}

	public static byte[] getEndSoftClipBaseQualities(SAMRecord record) {
		byte[] seq = record.getBaseQualities();
		if (seq == null)  return null;
		if (seq == SAMRecord.NULL_QUALS) return null;
		seq = Arrays.copyOfRange(seq, record.getReadLength() - getEndSoftClipLength(record), record.getReadLength());
		return seq;
	}

	public static SAMRecord ensureNmTag(ReferenceSequenceFile ref, SAMRecord record) {
		if (record == null || record.getIntegerAttribute(SAMTag.NM.name()) != null) {
			return record;
		}
		recalculateTagNm(ref, record);
		return record;
	}

	public static int recalculateTagNm(ReferenceSequenceFile ref, SAMRecord record) {
		if (record == null) return 0;
		if (record.getReadUnmappedFlag() || record.getReadBases() == null || record.getReadBases() == SAMRecord.NULL_SEQUENCE) {
			record.setAttribute(SAMTag.NM.name(), null);
			return 0;
		} else {
			byte[] refSeq = ref
					.getSubsequenceAt(record.getReferenceName(), record.getAlignmentStart(), record.getAlignmentEnd())
					.getBases();
			final int actualNucleotideDiffs = SequenceUtil.calculateSamNmTag(record, refSeq,
					record.getAlignmentStart() - 1);
			record.setAttribute(SAMTag.NM.name(), actualNucleotideDiffs);
			return actualNucleotideDiffs;
		}
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
		if (record1.getReadUnmappedFlag() || record2.getReadUnmappedFlag()) return false;
		return isDovetailing(record1.getReferenceIndex(), record1.getAlignmentStart(),
				record1.getReadNegativeStrandFlag(), record1.getCigar(), record2.getReferenceIndex(),
				record2.getAlignmentStart(), record2.getReadNegativeStrandFlag(), record2.getCigar(),
				expectedOrientation, margin);
	}

	public static boolean isDovetailing(SAMRecord record, PairOrientation expectedOrientation, int margin) {
		if (record.getReadUnmappedFlag()) return false;
		if (!record.getReadPairedFlag()) return false;
		if (record.getMateUnmappedFlag()) return false;
		Cigar cigar2 = SAMRecordUtil.getCachedMateCigar(record);
		return isDovetailing(record.getReferenceIndex(), record.getAlignmentStart(), record.getReadNegativeStrandFlag(),
				record.getCigar(), record.getMateReferenceIndex(), record.getMateAlignmentStart(),
				record.getMateNegativeStrandFlag(), cigar2, expectedOrientation, margin);
	}

	private static boolean isDovetailing(int reference1, int start1, boolean isNegativeStrand1, Cigar cigar1,
			int reference2, int start2, boolean isNegativeStrand2, Cigar cigar2, PairOrientation expectedOrientation,
			int margin) {
		if (expectedOrientation != PairOrientation.FR) throw new RuntimeException("NYI");
		if (reference1 != reference2) return false;
		if (Math.abs(start1 - start2) > margin) return false;
		if (isNegativeStrand1 == isNegativeStrand2) return false; // FR
		if (cigar2 != null) {
			int end1 = start1 + CigarUtil.referenceLength(cigar1.getCigarElements()) - 1;
			int end2 = start2 + CigarUtil.referenceLength(cigar2.getCigarElements()) - 1;
			if (Math.abs(end1 - end2) > margin) return false;
			if (!IntervalUtil.overlapsClosed(start1, end1, start2, end2)) return false;
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
		if (unexpectedClipLength > margin) return false;
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
		Cigar mc = SAMRecordUtil.getCachedMateCigar(record);
		if (record.getReadNegativeStrandFlag()) {
			// <--record
			int r1end = record.getUnclippedEnd();
			int r2start = mc == null ?
					// if we don't have a mate cigar we'll just assume that there are is no clipping in the alignment
					record.getMateAlignmentStart() :
					SAMUtils.getUnclippedStart(record.getMateAlignmentStart(), mc);
			return r1end - r2start + 1;
		} else {
			int r1start = record.getUnclippedStart();
			int r2end = mc == null ?
					// no MC tag: we have to assume no clipping and the reads are the same length
					record.getMateAlignmentStart() + CigarUtil.readLength(record.getCigar().getCigarElements()) + CigarUtil.countBases(record.getCigar(), CigarOperator.HARD_CLIP) - 1:
					SAMUtils.getUnclippedEnd(record.getMateAlignmentStart() + mc.getReferenceLength() - 1, mc);
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
		if (read.getReadBases() == null || read.getReadBases().length == 0) return 4;
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
	 * Recalculates the FI tag.
	 */
	public static final Set<String> TEMPLATE_TAGS = ImmutableSet.of(SAMTag.CC.name(), SAMTag.CP.name(), SAMTag.FI.name(), SAMTag.HI.name(),
			SAMTag.IH.name(), SAMTag.Q2.name(), SAMTag.R2.name(), SAMTag.MC.name(), SAMTag.MQ.name(), SAMTag.SA.name(), SAMTag.TC.name());
	public static List<List<SAMRecord>> calculateTagFI(List<List<SAMRecord>> templateBySegments) {
		// FI
		for (int i = 0; i < templateBySegments.size(); i++) {
			for (SAMRecord r : templateBySegments.get(i)) {
				r.setAttribute(SAMTag.FI.name(), i);
			}
		}
		return templateBySegments;
	}

	/**
	 * Recalculates the TC tag.
	 */
	public static List<List<SAMRecord>> calculateTagTC(List<List<SAMRecord>> templateBySegments) {
		// TC
		for (List<SAMRecord> segment : templateBySegments) {
			for (SAMRecord r : segment) {
				r.setAttribute(SAMTag.TC.name(), templateBySegments.size());
			}
		}
		return templateBySegments;
	}
	public static List<List<SAMRecord>> calculateTagR2(List<List<SAMRecord>> templateBySegments) {
		for (int i = 0; i < templateBySegments.size(); i++) {
			byte[] br2 = getConsensusSequence(templateBySegments.get((i + 1) % templateBySegments.size()));
			String r2 = br2 != null && br2 != SAMRecord.NULL_SEQUENCE && templateBySegments.size() > 1 ? StringUtil.bytesToString(br2) : null;
			for (SAMRecord r : templateBySegments.get(i)) {
				r.setAttribute(SAMTag.R2.name(), r2);
			}
		}
		return templateBySegments;
	}
	public static List<List<SAMRecord>> calculateTagQ2(List<List<SAMRecord>> templateBySegments) {
		for (int i = 0; i < templateBySegments.size(); i++) {
			byte[] bq2 = getConsensusBaseQualities(templateBySegments.get((i + 1) % templateBySegments.size()));
			String q2 = bq2 != null && bq2 != SAMRecord.NULL_QUALS && templateBySegments.size() > 1 ? SAMUtils.phredToFastq(bq2) : null;
			for (SAMRecord r : templateBySegments.get(i)) {
				r.setAttribute(SAMTag.Q2.name(), q2);
			}
		}
		return templateBySegments;
	}
	/**
	 * Ensures the duplicate flag is consistent.
	 *
	 * If any record in the template is flagged as a duplicate then all all.
	 *
	 * This fixes the issue in which duplicate marking software does not mark supplementary alignments
	 * as duplicates.
	 *
	 * @param templateRecords all records for a given template
	 * @return records
	 */
	public static List<SAMRecord> ensureConsistentDuplicateFlag(List<SAMRecord> templateRecords) {
		boolean isDuplicate = false;
		for (SAMRecord r : templateRecords) {
			if (r.getDuplicateReadFlag()) {
				isDuplicate = true;
				break;
			}
		}
		if (isDuplicate) {
			for (SAMRecord r : templateRecords) {
				r.setDuplicateReadFlag(true);
			}
		}
		return templateRecords;
	}
	private static final Comparator<SAMRecord> ByReadLength = Comparator.comparingInt(r -> r.getReadLength());
	public static List<SAMRecord> addMissingHardClipping(List<SAMRecord> segmentRecords) {
		if (segmentRecords.isEmpty()) return segmentRecords;
		SAMRecord longest = segmentRecords.stream().max(ByReadLength).get();
		String longseq = longest.getReadString();
		if (longest.getReadNegativeStrandFlag()) {
			longseq = SequenceUtil.reverseComplement(longseq).toUpperCase();
		}
		for (SAMRecord r : segmentRecords) {
			if (r.getReadUnmappedFlag()) continue;
			if (r.getReadLength() == longseq.length()) continue;
			List<CigarElement> cigar = r.getCigar().getCigarElements();
			int cigarLengthWithHardClips = r.getCigar().getReadLength() + cigar.stream()
					.filter(ce -> ce.getOperator() == CigarOperator.HARD_CLIP)
					.mapToInt(ce -> ce.getLength())
					.sum();
			if (cigarLengthWithHardClips == longseq.length()) continue;
			String seq = r.getReadString().toUpperCase();
			String fullseq = r.getReadNegativeStrandFlag() ? SequenceUtil.reverseComplement(longseq) : longseq;
			boolean isFixed = false;
			for (int i = 0 ; i <= fullseq.length() - seq.length(); i++) {
				String matchseq = fullseq.substring(i, i + seq.length());
				if (seq.equals(matchseq)) {
					List<CigarElement> newCigar = cigar.stream().filter(ce -> ce.getOperator() != CigarOperator.HARD_CLIP).collect(Collectors.toList());
					int leftPad = i;
					int rightPad = fullseq.length() - seq.length() - i;
					if (leftPad > 0) {
						newCigar.add(0, new CigarElement(leftPad, CigarOperator.HARD_CLIP));
					}
					if (rightPad > 0) {
						newCigar.add(new CigarElement(rightPad, CigarOperator.HARD_CLIP));
					}
					isFixed = true;
					r.setCigar(new Cigar(newCigar));
					break;
				}
			}
			if (!isFixed) {
				if (!MessageThrottler.Current.shouldSupress(log, "Add missing hard clipping")) {
					log.warn(String.format("Unable to restore hard clipping of truncated read %s. No sequence match to full length alignment record.", segmentRecords.get(0).getReadName()));
				}
			}
		}
		return segmentRecords;
	}

	/**
	 * Orders the records such that the primary record for a split read alignment is first
	 */
	private static Ordering<SAMRecord> ByBestPrimarySplitCandidate = new Ordering<SAMRecord>() {
		public int compare(SAMRecord arg1, SAMRecord arg2) {
			return ComparisonChain.start()
				.compareFalseFirst(arg1.getReadUnmappedFlag(), arg2.getReadUnmappedFlag())
				// already flagged as supp is bad  
				.compareFalseFirst(arg1.getSupplementaryAlignmentFlag(), arg2.getSupplementaryAlignmentFlag())
				// flagged as secondary is bad due to legacy treatment of secondary alignments as supplementary (eg bwa mem -M) 
				.compareFalseFirst(arg1.isSecondaryAlignment(), arg2.isSecondaryAlignment())
				// aligning to likely alt contigs is bad (https://github.com/lh3/bwa/issues/282)
				.compareFalseFirst(isLikelyAltContig(arg1), isLikelyAltContig(arg2))
				// longer read is better
				.compare(arg2.getReadLength(), arg1.getReadLength())
				// high MAPQ is better
				.compare(arg2.getMappingQuality(), arg1.getMappingQuality())
				// the record with the shorter soft clip is a better candidate
				.compare(
					SAMRecordUtil.getStartClipLength(arg1) + SAMRecordUtil.getEndClipLength(arg1),
					SAMRecordUtil.getStartClipLength(arg2) + SAMRecordUtil.getEndClipLength(arg2))
				.compare(
					Math.max(SAMRecordUtil.getStartClipLength(arg1), SAMRecordUtil.getEndClipLength(arg1)),
					Math.max(SAMRecordUtil.getStartClipLength(arg2), SAMRecordUtil.getEndClipLength(arg2)))
				// in proper pair is better
				.compareTrueFirst(arg1.getReadPairedFlag() && arg1.getProperPairFlag(), arg2.getReadPairedFlag() && arg2.getProperPairFlag())
				// lower edit distance is better
				.compare(getNM(arg1, 0), getNM(arg2, 0))
				// Other options are:
				// - record is flagged as the mate of a read pair (strong support for that record to be the primary)
				// Force a stable record sort order
				.compare(arg1.getReferenceIndex(), arg2.getReferenceIndex())
				.compare(arg1.getAlignmentStart(), arg2.getAlignmentStart())
				.compare(arg1.getCigarString(), arg2.getCigarString())
				.result();
		}
	};
	private static boolean isLikelyAltContig(SAMRecord r) {
		if (r.getReadUnmappedFlag()) return false;
		// hg19
		if (r.getReferenceName().contains("_random")) return true;
		if (r.getReferenceName().contains("Un_gl")) return true;
		// hg38
		if (r.getReferenceName().contains("decoy")) return true;
		if (r.getReferenceName().contains("HLA-")) return true;
		if (r.getReferenceName().contains("_alt")) return true;
		return false;
	}

	private static int getNM(SAMRecord r, int defaultValue) {
		Integer value = r.getIntegerAttribute(SAMTag.NM.name());
		return value == null ? defaultValue : value;
	}
	public static void recalculateSupplementaryFromSA(List<SAMRecord> segments) {
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
	public static boolean matchReadPairPrimaryAlignments(List<List<SAMRecord>> segments) {
		if (segments.size() != 2) return false;
		List<SAMRecord> segment1 = segments.get(0);
		List<SAMRecord> segment2 = segments.get(1);
		if (segment1.size() == 0) return false;
		if (segment2.size() == 0) return false;
		List<SAMRecord> primary1 = segment1.stream().filter(r -> !r.isSecondaryOrSupplementary()).collect(Collectors.toList());
		List<SAMRecord> primary2 = segment2.stream().filter(r -> !r.isSecondaryOrSupplementary()).collect(Collectors.toList());
		if (primary1.size() == 1 && primary2.size() == 1) {
			SAMRecord r1 = primary1.get(0);
			SAMRecord r2 = primary2.get(0);
			SamPairUtil.setMateInfo(r1, r2, true);
			r1.removeTransientAttribute(SAMTag.MC.name());
			r2.removeTransientAttribute(SAMTag.MC.name());
		}
		return true;
	}

	public static void fixMates(List<List<SAMRecord>> segments, boolean mateCigar, boolean mateQuality) {
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
				r.removeTransientAttribute(SAMTag.MC.name());
				if (r.getSupplementaryAlignmentFlag()) {
					SAMRecord primary = getPrimarySplitAlignmentFor(r, currentSegment);
					SAMRecord mate = bestMateFor(primary, nextSegment);
					if (mate != null) {
						SamPairUtil.setMateInformationOnSupplementalAlignment(r, mate, mateCigar);
						if (mateQuality) {
							r.setAttribute(SAMTag.MQ.name(), mate.getMappingQuality());
						}
					} else {
						clearMateInformation(r, true);
					}
				} else {
					SAMRecord mate = bestMateFor(r, nextSegment);
					if (mate != null) {
						SamPairUtil.setMateInfo(r, mate, mateCigar);
						if (mateQuality) {
							r.setAttribute(SAMTag.MQ.name(), mate.getMappingQuality());
						}
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
	public static void clearMateInformation(SAMRecord r, boolean fullClear) {
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
			r.removeTransientAttribute(SAMTag.MC.name());
		}
	}

	private static SAMRecord bestMateFor(SAMRecord r, List<SAMRecord> mates) {
		// Mates always point to supp alignments
		List<SAMRecord> candidates = mates.stream()
				.filter(m -> !m.getSupplementaryAlignmentFlag())
				.collect(Collectors.toList());
		// see if we already point to one of the mates (ie: don't touch existing records)
		Optional<SAMRecord> best = candidates.stream()
				.filter(m -> m.getReadUnmappedFlag() == r.getMateUnmappedFlag())
				.filter(m -> m.getReferenceIndex() == r.getMateReferenceIndex())
				.filter(m -> m.getAlignmentStart() == r.getMateAlignmentStart())
				.filter(m -> m.getReadNegativeStrandFlag() == r.getMateNegativeStrandFlag())
				.filter(m -> m.getCigarString().equals(r.getStringAttribute("MC")) || r.getStringAttribute("MC")== null)
				.filter(m -> m.isSecondaryAlignment() == r.isSecondaryAlignment())
				.findFirst();
		if (best.isPresent()) return best.get();
		if (!r.getSupplementaryAlignmentFlag()) {
			// see if one of the mates already points to us
			best = candidates.stream()
					.filter(m -> r.getReadUnmappedFlag() == m.getMateUnmappedFlag())
					.filter(m -> r.getReferenceIndex() == m.getMateReferenceIndex())
					.filter(m -> r.getAlignmentStart() == m.getMateAlignmentStart())
					.filter(m -> r.getReadNegativeStrandFlag() == m.getMateNegativeStrandFlag())
					.findFirst();
			if (best.isPresent()) return best.get();
		}
		// match primary with primary and secondary with secondary
		best = candidates.stream()
				.filter(m -> m.isSecondaryAlignment() == r.isSecondaryAlignment())
				// prefer mapped reads
				.sorted(Comparator.comparing(SAMRecord::getReadUnmappedFlag)).findFirst();
		if (best.isPresent()) return best.get();
		// just get anything that's not a supplementary alignment
		best = mates.stream()
				.filter(m -> m.getSupplementaryAlignmentFlag() == r.getSupplementaryAlignmentFlag())
				.findFirst();
		if (best.isPresent()) return best.get();
		// we give up - anything will do at this point
		best = mates.stream().findFirst();
		return best.orElse(null);
	}

	/**
	 * Groups records according to their segment index.
	 *
	 * List lengths is 1 for unpaired reads, 2 for paired reads.
	 */
	public static final List<List<SAMRecord>> templateBySegment(List<SAMRecord> records) {
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
	// TODO support secondary alignments via CC, CP, HI, IH
	public static List<SAMRecord> reinterpretAsSplitReadAlignment(List<SAMRecord> segmentRecords, int overlapThreshold) {
		if (segmentRecords == null || segmentRecords.isEmpty()) {
			return Collections.emptyList();
		}
		try {
			segmentRecords.sort(ByBestPrimarySplitCandidate);
		} catch (UnsupportedOperationException e) {
			segmentRecords = Lists.newArrayList(segmentRecords);
			segmentRecords.sort(ByBestPrimarySplitCandidate);
		}
		SAMRecord primary = segmentRecords.get(0);
		// TODO: actually support secondary alignments
		for (SAMRecord r : segmentRecords) {
			r.setSecondaryAlignment(false);
			r.setSupplementaryAlignmentFlag(r != primary);
		}
		// Remove records that overlap or contain 'better' alignments
		for (int i = 1; i < segmentRecords.size(); i++) {
			SAMRecord r = segmentRecords.get(i);
			int start = getFirstAlignedBaseReadOffset(r);
			int end = getLastAlignedBaseReadOffset(r);
			int totalOverlap = 0;
			for (int j = 0; j < i; j++) {
				SAMRecord prev = segmentRecords.get(j);
				int prevStart = getFirstAlignedBaseReadOffset(prev);
				int prevEnd = getLastAlignedBaseReadOffset(prev);
				int minLength = Math.min(prevEnd - prevStart, end - start) + 1;
				int overlap = IntervalUtil.overlapsWidthClosed(start, end, prevStart, prevEnd);
				totalOverlap += overlap;
				if (overlap == minLength || totalOverlap > overlapThreshold) {
					segmentRecords.remove(i);
					i--;
					break;
				}
			}
		}
		if (segmentRecords.size() == 1) {
			segmentRecords.get(0).setAttribute(SAMTag.SA.name(), null);
		} else {
			// primary is at start of list
			List<String> saString = new ArrayList<>(segmentRecords.size());
			for (SAMRecord r : segmentRecords) {
				StringBuilder sb = new StringBuilder();
				for (SAMRecord r2 : segmentRecords) {
					if (r != r2) {
						sb.append(new ChimericAlignment(r2).toString());
						sb.append(';');
					}
				}
				// strip trailing semicolon
				sb.setLength(sb.length() - 1);
				r.setAttribute(SAMTag.SA.name(), sb.toString());
			}
		}
		return segmentRecords;
	}

	/**
	 * replaced by reinterpretAsSplitReadAlignment()
	 */
	@Deprecated
	private static void calculateNonoverlappingSplitReadAlignmentSATags(List<SAMRecord> list, boolean recalculateSA) {
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
			if (!recalculateSA && r.getAttribute(SAMTag.SA.name()) != null) {
				continue;
			}
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
		warnIfInvalidSA(list);
	}

	private static void warnIfInvalidSA(List<SAMRecord> list) {
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
					r.setTransientAttribute("HC", null);
					if (qual != null && qual != SAMRecord.NULL_QUALS) {
						byte[] newqual = Arrays.copyOf(qual, qual.length);
						if (r.getReadNegativeStrandFlag()) {
							ArrayUtils.reverse(newqual);
						}
						r.setBaseQualities(newqual);
					}
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
		r.setTransientAttribute("HC", true);
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
		return getFirstAlignedBaseReadOffset(r.getCigar(), r.getReadNegativeStrandFlag());

	}
	public static final int getFirstAlignedBaseReadOffset(Cigar c, boolean onNegativeStrand) {
		if (c == null || c.getCigarElements().size() == 0)
			return -1;
		if (onNegativeStrand) {
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
		return getLastAlignedBaseReadOffset(r.getCigar(), r.getReadNegativeStrandFlag());
	}
	public static final int getLastAlignedBaseReadOffset(Cigar c, boolean onNegativeStrand) {
		if (c == null || c.getCigarElements().size() == 0)
			return -1;
		int length = c.getReadLength();
		if (c.getFirstCigarElement().getOperator() == CigarOperator.HARD_CLIP) {
			length += c.getFirstCigarElement().getLength();
		}
		if (c.getLastCigarElement().getOperator() == CigarOperator.HARD_CLIP) {
			length += c.getLastCigarElement().getLength();
		}
		if (onNegativeStrand) {
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
		adjustAlignmentBounds(read, startunclip, endunclip);
	}
	/**
	 * Adjusts the start/end alignment bounds of the given alignment. Expanded alignments are converted
	 * @param r
	 * @param expandStartBy
	 * @param expandEndBy
	 * @return record passed in
	 */
	public static SAMRecord adjustAlignmentBounds(SAMRecord r, int expandStartBy, int expandEndBy) {
		Cigar newCigar = softClipToAligned(r.getCigar(), expandStartBy, 0);
		r.setAlignmentStart(r.getAlignmentStart() + r.getCigar().getReferenceLength() - newCigar.getReferenceLength());
		newCigar = softClipToAligned(newCigar, 0, expandEndBy);
		r.setCigar(newCigar);
		return r;
	}
	private static Cigar softClipToAligned(Cigar cigar, int startunclip, int endunclip) {
		if (startunclip == 0 && endunclip == 0) {
			return cigar;
		}
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
		} else if (startunclip < 0) {
			int i = 0;
			int alignedBasesLeftToRemove = -startunclip;
			while (!list.get(i).getOperator().consumesReferenceBases()) {
				i++;
			}
			while (alignedBasesLeftToRemove > 0 && i >= 0 && i < list.size()) {
				CigarElement ce = list.get(i);
				if (ce.getOperator().consumesReadBases()) {
					if (ce.getLength() >= alignedBasesLeftToRemove) {
						list.remove(i);
						if (ce.getLength() != alignedBasesLeftToRemove) {
							list.add(i, new CigarElement(ce.getLength() - alignedBasesLeftToRemove, ce.getOperator()));
						}
						list.add(i, new CigarElement(alignedBasesLeftToRemove, CigarOperator.SOFT_CLIP));
						alignedBasesLeftToRemove = 0;
					} else {
						list.set(i, new CigarElement(ce.getLength(), CigarOperator.SOFT_CLIP));
						alignedBasesLeftToRemove -= ce.getLength();
						i++;
					}
				} else {
					// D and N operators do not consume read bases - throw them out
					list.remove(i);
				}
			}
			CigarUtil.clean(list, false);
		}
		if (endunclip > 0) {
			int i = list.size() - 1;
			while (i >= 0 && list.get(i).getOperator() != CigarOperator.SOFT_CLIP) {
				i--;
			}
			if (i < 0 || list.get(i).getLength() < endunclip) {
				throw new IllegalArgumentException(
						String.format("Unable to unclip %d clipped bases from end of %s", endunclip, cigar));
			}
			list.set(i, new CigarElement(list.get(i).getLength() - endunclip, CigarOperator.SOFT_CLIP));
			list.add(i - 1, new CigarElement(endunclip, CigarOperator.MATCH_OR_MISMATCH));
			CigarUtil.clean(list, false);
		} else if (endunclip < 0) {
			int i = list.size() - 1;
			while (i >= 0 && !list.get(i).getOperator().consumesReferenceBases()) {
				i--;
			}
			int alignedBasesLeftToRemove = -endunclip;
			while (alignedBasesLeftToRemove > 0 && i >= 0 && i < list.size()) {
				CigarElement ce = list.get(i);
				if (ce.getOperator().consumesReadBases()) {
					if (ce.getLength() >= alignedBasesLeftToRemove) {
						list.remove(i);
						list.add(i, new CigarElement(alignedBasesLeftToRemove, CigarOperator.SOFT_CLIP));
						if (ce.getLength() != alignedBasesLeftToRemove) {
							list.add(i, new CigarElement(ce.getLength() - alignedBasesLeftToRemove, ce.getOperator()));
						}
						alignedBasesLeftToRemove = 0;
					} else {
						alignedBasesLeftToRemove -= ce.getLength();
						list.set(i, new CigarElement(ce.getLength(), CigarOperator.SOFT_CLIP));
					}
				} else {
					// D and N operators do not consume read bases - throw them out
					list.remove(i);
				}
				i--;
			}
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
	/**
	 * Returns the number of bases common between the given alignments
	 * @param r1 record to compare
	 * @param r2 record to compare
	 * @return number of matching reference-aligned bases.
	 */
	public static int overlappingBases(SAMRecord r1, SAMRecord r2) {
		if (!overlap(r1, r2)) {
			return 0;
		}
		return overlappingBases(
				r1.getReferenceIndex(), r1.getAlignmentStart(), r1.getReadNegativeStrandFlag(), r1.getCigar(),
				r2.getReferenceIndex(), r2.getAlignmentStart(), r2.getReadNegativeStrandFlag(), r2.getCigar());
	}
	public static int overlappingBases(
			int referenceIndex1, int alignmentStart1, boolean negativeStrand1, Cigar cigar1,
			int referenceIndex2, int alignmentStart2, boolean negativeStrand2, Cigar cigar2) {
		if (referenceIndex1 != referenceIndex2 || negativeStrand1 != negativeStrand2) {
			return 0;
		}
		ImmutableRangeSet<Integer> rs1 = getMappedBases(alignmentStart1, cigar1);
		ImmutableRangeSet<Integer> rs2 = getMappedBases(alignmentStart2, cigar2);
		ImmutableRangeSet<Integer> overlapSet = rs1.intersection(rs2);
		int overlapCount = 0;
		for (Range<Integer> r : overlapSet.asRanges()) {
			overlapCount += r.upperEndpoint() - r.lowerEndpoint();
		}
		return overlapCount;
	}
	private static ImmutableRangeSet<Integer> getMappedBases(int start, Cigar cigar) {
		Builder<Integer> builder = ImmutableRangeSet.builder();
		int position = start;
		for (CigarElement op : cigar.getCigarElements()) {
			if (op.getOperator().consumesReferenceBases()) {
				if (op.getOperator().consumesReadBases()) {
					builder.add(Range.closedOpen(position, position + op.getLength()));
				}
				position += op.getLength();
			}
		}
		return builder.build();
	}
	/**
	 * Determines whether any read alignments overlaps the original alignment (if any).
	 * @param r
	 * @return
	 */
	public static boolean overlapsOriginalAlignment(SAMRecord r) {
		String oa = r.getStringAttribute("OA");
		if (StringUtil.isBlank(oa)) return false;
		SAMSequenceDictionary dict = r.getHeader().getSequenceDictionary();
		ChimericAlignment originalAlignment = new ChimericAlignment(oa);
		ChimericAlignment thisAlignment = new ChimericAlignment(r);
		int overlapCount = originalAlignment.overlappingBases(dict, thisAlignment);
		for (ChimericAlignment ca : ChimericAlignment.getChimericAlignments(r)) {
			overlapCount += originalAlignment.overlappingBases(dict, ca);
		}
		return overlapCount > 0;
	}
	public static int getReadLengthIncludingHardClipping(SAMRecord r) {
		final int lengthWithHardClipping = r.getReadLength() + r.getCigar().getCigarElements().stream()
			.filter(ce -> ce.getOperator() == CigarOperator.HARD_CLIP)
			.mapToInt(ce -> ce.getLength()).sum();
		return lengthWithHardClipping;
	}
	public static Cigar getCachedMateCigar(SAMRecord r) {
		Object cached = r.getTransientAttribute(SAMTag.MC.name());
		if (cached != null) {
			return (Cigar)cached;
		} else {
			String s = r.getStringAttribute(SAMTag.MC.name());
			if (s != null) {
				Cigar cigar = TextCigarCodec.decode(s);
				r.setTransientAttribute(cigar, SAMTag.MC.name());
				return cigar;
			}
		}
		return null;
	}
	public static boolean forceValidContigBounds(SAMRecord r, SAMSequenceDictionary dict) {
		if (r.getReadUnmappedFlag()) return false;
		if (dict == null) return false;
		if (r.getReferenceIndex() == null) return false;
		if (dict.getSequence(r.getReferenceIndex()) == null) return false;
		int seqlen = Integer.MAX_VALUE;
		if (dict != null && r.getReferenceIndex() != null && dict.getSequence(r.getReferenceIndex()) != null) {
			seqlen = dict.getSequence(r.getReferenceIndex()).getSequenceLength();
		}
		if (r.getAlignmentStart() > 0 && r.getAlignmentEnd() <= seqlen) return false;
		// trim the start
		int endPosition = r.getAlignmentEnd();
		List<CigarElement> cigar = forceStartOnContig(r.getAlignmentStart(), r.getCigar().getCigarElements());
		r.setCigar(new Cigar(cigar));
		int newEndPosition = r.getAlignmentEnd();
		// adjust start so we're still ending at the same position
		r.setAlignmentStart(r.getAlignmentStart() + endPosition - newEndPosition);
		// trim the end by flipping the cigar and reusing the start trimming code
		cigar = Lists.reverse(cigar);
		cigar = forceStartOnContig(seqlen - r.getAlignmentEnd() + 1, cigar);
		cigar = Lists.reverse(cigar);
		r.setCigar(new Cigar(CigarUtil.clean(cigar, false)));
		return false;
	}

	/**
	 * Adjusts the CIGAR so that the startPosition is at least 1
	 */
	private static List<CigarElement> forceStartOnContig(int startPosition, List<CigarElement> cigar) {
		List<CigarElement> out = new ArrayList<>();
		int positionsLeftToUnmap = 1 - startPosition;
		for (CigarElement ce : cigar) {
			if (positionsLeftToUnmap <= 0) {
				// pass through the rest of the cigar untouched
				out.add(ce);
			} else if (ce.getOperator().isClipping()) {
				out.add(ce);
			} else if (ce.getOperator().consumesReadBases() && ce.getOperator().consumesReferenceBases()) {
				int basesToConvertToSoftClip = Math.min(ce.getLength(), positionsLeftToUnmap);
				int basesToRetain = ce.getLength() - basesToConvertToSoftClip;
				if (basesToConvertToSoftClip > 0) {
					out.add(new CigarElement(basesToConvertToSoftClip, CigarOperator.SOFT_CLIP));
				}
				if (basesToRetain > 0) {
					out.add(new CigarElement(basesToRetain, ce.getOperator()));
				}
				positionsLeftToUnmap -= basesToConvertToSoftClip;
			} else if (ce.getOperator().consumesReadBases() && !ce.getOperator().consumesReferenceBases()) {
				out.add(new CigarElement(ce.getLength(), CigarOperator.SOFT_CLIP));
			} else if (!ce.getOperator().consumesReadBases() && ce.getOperator().consumesReferenceBases()) {
				positionsLeftToUnmap -= ce.getLength();
			} else {
				// doesn't consume read or reference bases - strip
			}
		}
		return out;
	}

	public static SAMRecord createSAMRecord(SAMFileHeader header, FastqRecord fq, boolean reverseComp) {
		SAMRecord r = new SAMRecord(header);
		byte[] seq = fq.getReadBases();
		byte[] qual = fq.getBaseQualities();
		if (reverseComp) {
			SequenceUtil.reverseComplement(seq);
			ArrayUtils.reverse(qual);
		}
		r.setReadName(fq.getReadName());
		r.setReadBases(seq);
		r.setBaseQualities(qual);
		return r;
	}

	/**
	 * Determines whether the alignments overlap on the same strand
	 * @param r1 first record
	 * @param r2 second record
	 * @return true if the alignments overlap, false otherwise
	 */
	public static boolean alignmentOverlaps(SAMRecord r1, SAMRecord r2) {
		return !r1.getReadUnmappedFlag() && !r2.getReadUnmappedFlag() && r1.getReadNegativeStrandFlag() == r2.getReadNegativeStrandFlag() && r1.overlaps(r2);
	}

	private static final Ordering<SAMRecord> ByFirstAlignedBaseReadOffset = new Ordering<SAMRecord>() {
		public int compare(SAMRecord arg1, SAMRecord arg2) {
			return Ints.compare(getFirstAlignedBaseReadOffset(arg1), getFirstAlignedBaseReadOffset(arg2));
		}
	};

	public static boolean clipTerminalIndelCigars(SAMRecord r) {
		if (r.getReadUnmappedFlag()) return false;
		Cigar c = r.getCigar();
		if (!c.containsOperator(CigarOperator.INSERTION) && !c.containsOperator(CigarOperator.DELETION)) return false;
		List<CigarElement> cigar = Lists.newArrayList(r.getCigar().getCigarElements());
		Pair<Integer, Integer> endBasesAdjusted = clipIndelEndCigar(cigar);
		Pair<Integer, Integer> startBasesAdjusted = clipIndelStartCigar(cigar);
		if (endBasesAdjusted.getFirst() + endBasesAdjusted.getSecond() + startBasesAdjusted.getFirst() + startBasesAdjusted.getSecond() == 0) return false;
		int alignmentOffset = startBasesAdjusted.getSecond();
		CigarUtil.clean(cigar, false);
		r.setCigar(new Cigar(cigar)); // NB: https://github.com/samtools/htsjdk/pull/1538
		r.setAlignmentStart(r.getAlignmentStart() + alignmentOffset);
		boolean hasAlignedBase = cigar.stream().anyMatch(ce -> ce.getOperator().consumesReferenceBases());
		r.setReadUnmappedFlag(!hasAlignedBase);
		return true;
	}

	/**
	 *
	 * @return bases removed from (read, reference)
	 */
	private static Pair<Integer, Integer> clipIndelEndCigar(List<CigarElement> cigar) {
		int insBases = 0;
		int delBases = 0;
		for (int i = cigar.size() - 1; i >= 0; i--) {
			CigarElement ce = cigar.get(i);
			if (ce.getOperator().isIndel()) {
				if (ce.getOperator() == CigarOperator.INSERTION) {
					cigar.set(i, new CigarElement(ce.getLength(), CigarOperator.SOFT_CLIP));
					insBases += ce.getLength();
				} else {
					cigar.remove(i);
					delBases += ce.getLength();
				}
			} else if (!ce.getOperator().isClipping()) {
				break;
			}
		}
		return Pair.create(insBases, delBases);
	}
	private static Pair<Integer, Integer> clipIndelStartCigar(List<CigarElement> cigar) {
		int insBases = 0;
		int delBases = 0;
		for (int i = 0; i < cigar.size(); i++) {
			CigarElement ce = cigar.get(i);
			if (ce.getOperator().isIndel()) {
				if (ce.getOperator() == CigarOperator.INSERTION) {
					cigar.set(i, new CigarElement(ce.getLength(), CigarOperator.SOFT_CLIP));
					insBases += ce.getLength();
				} else {
					cigar.remove(i);
					i--;
					delBases += ce.getLength();
				}
			} else if (!ce.getOperator().isClipping()) {
				break;
			}
		}
		return Pair.create(insBases, delBases);
	}
	/**
	 * Input records must be grouped by read name
	 */
	public static Iterator<List<SAMRecord>> groupedByReadName(Iterator<SAMRecord> it) {
		return new GroupingIterator(it, Ordering.natural().onResultOf((SAMRecord r) -> r.getReadName()));
	}

    public static boolean hasAlignment(SAMRecord record) {
		if (record.getReadUnmappedFlag()) return false;
		for (CigarElement ce : record.getCigar().getCigarElements()) {
			switch (ce.getOperator()) {
				case EQ:
				case M:
				case X:
					return true;
			}
		}
		return false;
    }
}

