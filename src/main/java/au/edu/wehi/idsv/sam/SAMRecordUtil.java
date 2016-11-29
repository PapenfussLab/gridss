/**
 * 
 */
package au.edu.wehi.idsv.sam;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.stream.Collectors;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.Pair;

import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;
import com.google.common.collect.Sets;
import com.google.common.primitives.Booleans;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.alignment.AlignerFactory;
import au.edu.wehi.idsv.alignment.Alignment;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.util.IntervalUtil;
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
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;

/**
 * @author Daniel Cameron
 *
 */
public class SAMRecordUtil {
	private static final Log log = Log.getInstance(SAMRecordUtil.class);
	private static final String SUFFIX_SEPERATOR = "#";
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
		if (cigar == null) return 0;
		return getStartClipLength(cigar.getCigarElements());
	}
	
	public static int getStartClipLength(List<CigarElement> elements) {
		if (elements == null) return 0;
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
	
	public static int getEndClipLength(SAMRecord aln) {
		Cigar cigar = aln.getCigar();
		if (cigar == null) return 0;
		return getEndClipLength(cigar.getCigarElements());
	}
	
	public static int getEndClipLength(List<CigarElement> elements) {
		if (elements == null) return 0;
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
		byte[] refSeq = ref.getSubsequenceAt(record.getReferenceName(), record.getAlignmentStart(), record.getAlignmentEnd()).getBases();
		final int actualNucleotideDiffs = SequenceUtil.calculateSamNmTag(record, refSeq, record.getAlignmentStart() - 1);
		record.setAttribute(SAMTag.NM.name(), actualNucleotideDiffs);
		return record;
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
		if (record1.getReadUnmappedFlag() || record2.getReadUnmappedFlag()) return false;
		return isDovetailing(
				record1.getReferenceIndex(), record1.getAlignmentStart(), record1.getReadNegativeStrandFlag(), record1.getCigar(),
				record2.getReferenceIndex(), record2.getAlignmentStart(), record2.getReadNegativeStrandFlag(), record2.getCigar(),
				expectedOrientation, margin);
	}
	public static boolean isDovetailing(SAMRecord record, PairOrientation expectedOrientation, int margin) {
		if (record.getReadUnmappedFlag()) return false;
		if (!record.getReadPairedFlag()) return false;
		if (record.getMateUnmappedFlag()) return false;
		Cigar cigar2 = null;
		String mc = record.getStringAttribute(SAMTag.MC.name());
		if (mc != null) {
			cigar2 = TextCigarCodec.decode(mc);
		}
		return isDovetailing(
				record.getReferenceIndex(), record.getAlignmentStart(), record.getReadNegativeStrandFlag(), record.getCigar(),
				record.getMateReferenceIndex(), record.getMateAlignmentStart(), record.getMateNegativeStrandFlag(), cigar2,
				expectedOrientation, margin);
	}
	private static boolean isDovetailing(
			int reference1,
			int start1,
			boolean isNegativeStrand1,
			Cigar cigar1,
			int reference2,
			int start2,
			boolean isNegativeStrand2,
			Cigar cigar2,
			PairOrientation expectedOrientation, 
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
		//    >>>SSS
		// SSS<<<
		// not like:
		// SSS>>>
		//    <<<SSS
		int unexpectedClipLength = isNegativeStrand1 ? getEndSoftClipLength(cigar1.getCigarElements()) : getStartSoftClipLength(cigar1.getCigarElements());
		if (cigar2 != null) {
			unexpectedClipLength += isNegativeStrand2 ? getEndSoftClipLength(cigar2.getCigarElements()) : getStartSoftClipLength(cigar2.getCigarElements());
		}
		if (unexpectedClipLength > margin) return false;
		return true;
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
		if (ass != null && ass.length > 0 && ref != null && ref.length > 0) {
			try {
				alignment = AlignerFactory.create().align_smith_waterman(ass, ref);        
			} catch (Exception e) {
				// swallow and log alignment error
				log.error(e, String.format("Error aligning %s to %s:%d-%d",
						read.getReadName(),
						refSeq.getSequenceName(), start, end));
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
	/**
	 * Computed tags requiring reads from the same template to be processed together
	 */
	public static final Set<SAMTag> TEMPLATE_TAGS = ImmutableSet.of(
			SAMTag.CC,
			SAMTag.CP,
			SAMTag.FI,
			SAMTag.HI,
			SAMTag.IH,
			SAMTag.Q2,
			SAMTag.R2,
			SAMTag.MC,
			SAMTag.MQ,
			SAMTag.SA,
			SAMTag.TC);
	/**
	 * Populates the missing template tags for the records originating from the same template 
	 * @param records Records from the same template.
	 * @param tags Tags to populate 
	 * @param restoreHardClips convert hard clips into soft clips if the sequence is available
	 * from a chimeric alignment of the same segmentt
	 * @param fixMates 
	 */
	public static void calculateTemplateTags(List<SAMRecord> records, Set<SAMTag> tags, boolean restoreHardClips, boolean fixMates) {
		List<List<SAMRecord>> segments = templateBySegment(records);
		// FI
		if (tags.contains(SAMTag.FI)) {
			for (int i = 0; i < segments.size(); i++) {
				for (SAMRecord r : segments.get(i)) {
					r.setAttribute(SAMTag.FI.name(), i);
				}
			}
		}
		// TC
		if (tags.contains(SAMTag.TC)) {
			for (SAMRecord r : records) {
				r.setAttribute(SAMTag.TC.name(), segments.size());
			}
		}
		// R2
		if (tags.contains(SAMTag.R2)) {
			for (int i = 0; i < segments.size(); i++) {
				byte[] br2 = getFullSequence(segments.get((i + 1) % segments.size()));
				String r2 = br2 != null && segments.size() > 1 ? StringUtil.bytesToString(br2) : null;
				for (SAMRecord r : segments.get(i)) {
					r.setAttribute(SAMTag.R2.name(), r2);
				}
			}
		}
		// Q2
		if (tags.contains(SAMTag.Q2)) {
			for (int i = 0; i < segments.size(); i++) {
				byte[] bq2 = getFullBaseQualities(segments.get((i + 1) % segments.size()));
				String q2 = bq2 != null && segments.size() > 1 ? SAMUtils.phredToFastq(bq2) : null;
				for (SAMRecord r : segments.get(i)) {
					r.setAttribute(SAMTag.Q2.name(), q2);
				}
			}
		}
		if (restoreHardClips) {
			for (int i = 0; i < segments.size(); i++) {
				softenHardClips(segments.get(i));
			}
		}
		// SA
		if (tags.contains(SAMTag.SA)) {
			for (int i = 0; i < segments.size(); i++) {
				calculateSATags(segments.get(i));
			}
		}
		if (Sets.intersection(tags, ImmutableSet.of(
				SAMTag.CC,
				SAMTag.CP,
				SAMTag.HI,
				SAMTag.IH)).size() > 0) {
			for (int i = 0; i < segments.size(); i++) {
				calculateMultimappingTags(tags, segments.get(i));
			}
		}
		if (fixMates || tags.contains(SAMTag.MC) || tags.contains(SAMTag.MQ)) {
			if (records.size() >= 2) {
				// hack to prevent SAMRecord getters from throwing an exception
				records.stream().forEach(r -> r.setReadPairedFlag(true));
				segments.get(0).stream().forEach(r -> r.setFirstOfPairFlag(true));
				segments.get(1).stream().forEach(r -> r.setSecondOfPairFlag(true));
				fixMates(segments, tags.contains(SAMTag.MC));
			} else {
				for (SAMRecord r : records) {
					clearMateInformation(r, true);
				}
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
					.filter(r -> r.getNotPrimaryAlignmentFlag() == rec.getNotPrimaryAlignmentFlag())
					.findFirst();
		}
		// grab anything that's not a supplementary
		if (!best.isPresent()) {
			best = options.stream()
					.filter(r -> !r.getSupplementaryAlignmentFlag())
					.findFirst();
		}
		// Just use ourself if there's nothing else out there but other supplementary alignments
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
			currentSegment.sort(Comparator.comparing(SAMRecord::getNotPrimaryAlignmentFlag).reversed());
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
	 * @param r read to clear mate information from
	 * @param fullClear strips all mate information from the read 
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
		// see if one of the mates already points to us
		Optional<SAMRecord> best = mates.stream()
				.filter(m -> r.getReadUnmappedFlag() == m.getMateUnmappedFlag())
				.filter(m -> r.getReferenceIndex() == m.getMateReferenceIndex())
				.filter(m -> r.getAlignmentStart() == m.getMateAlignmentStart())
				.filter(m -> r.getReadNegativeStrandFlag() == m.getMateNegativeStrandFlag())
				.findFirst();
		// see if we already point to one of the mates
		if (!best.isPresent()) {
			best = mates.stream()
					.filter(m -> m.getReadUnmappedFlag() == r.getMateUnmappedFlag())
					.filter(m -> m.getReferenceIndex() == r.getMateReferenceIndex())
					.filter(m -> m.getAlignmentStart() == r.getMateAlignmentStart())
					.filter(m -> m.getReadNegativeStrandFlag() == r.getMateNegativeStrandFlag())
					.findFirst();
		}
		// match primary with primary and secondary with secondary
		if (!best.isPresent()) {
			best = mates.stream()
					.filter(m -> m.getSupplementaryAlignmentFlag() == r.getSupplementaryAlignmentFlag())
					.filter(m -> m.getNotPrimaryAlignmentFlag() == r.getNotPrimaryAlignmentFlag())
					// prefer mapped reads
					.sorted(Comparator.comparing(SAMRecord::getReadUnmappedFlag))
					.findFirst();
		}
		// just get anything that's not a supplementary alignment
		if (!best.isPresent()) {
			best = mates.stream()
					.filter(m -> m.getSupplementaryAlignmentFlag() == r.getSupplementaryAlignmentFlag())
					.findFirst();
		}
		// we give up - anything will do at this point
		if (!best.isPresent()) {
			best = mates.stream()
					.findFirst();
		}
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
	public static void calculateMultimappingTags(Set<SAMTag> tags, List<SAMRecord> list) {
		// TODO: how does SA split read alignment interact with multimapping reads?
		list.sort(new SAMRecordCoordinateComparator());
		for (int i = 0; i < list.size(); i++) {
			SAMRecord r = list.get(i);
			if (tags.contains(SAMTag.IH)) {
				// TODO: how does this differ from NH?
				r.setAttribute(SAMTag.IH.name(), list.size());
			}
			if (tags.contains(SAMTag.HI)) {
				r.setAttribute(SAMTag.HI.name(), i);
			}
			if (tags.contains(SAMTag.CC)) {
				r.setAttribute(SAMTag.CC.name(), list.get((i + 1) % list.size()).getReferenceName());
			}
			if (tags.contains(SAMTag.CP)) {
				r.setAttribute(SAMTag.CP.name(), list.get((i + 1) % list.size()).getAlignmentStart());
			}
		}
	}

	private static void calculateSATags(List<SAMRecord> list) {
		// TODO use CC, CP, HI, IH tags if they are present
		// TODO break ties in an other other than just overwriting the previous alignment that
		// finishes at a given position
		SortedMap<Integer, SAMRecord> start = new TreeMap<Integer, SAMRecord>();
		SortedMap<Integer, SAMRecord> end = new TreeMap<Integer, SAMRecord>();
		for (SAMRecord r : list) {
			start.put(getFirstAlignedBaseReadOffset(r), r);
			end.put(getLastAlignedBaseReadOffset(r), r);
		}
		for (SAMRecord r : list) {
			if (r.getAttribute(SAMTag.SA.name()) != null) continue;
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
						// Conventionally, at a supplementary line, the first element points to the primary line.
						return Booleans.compare(arg0.getSupplementaryAlignmentFlag(), arg1.getSupplementaryAlignmentFlag());
					}});
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
			log.warn(String.format("SA tag of read %s refers to missing alignments %s", list.get(0).getReadName(), referencedReads));
		}
	}
	/**
	 * Converts any hard clips into soft clips using alternate alignments
	 * @param records
	 */
	public static final void softenHardClips(List<SAMRecord> records) {
		byte[] seq = getFullSequence(records);
		byte[] qual = getFullBaseQualities(records);
		if (seq == null) return;
		for (SAMRecord r : records) {
			Cigar c = r.getCigar();
			if (r.getReadUnmappedFlag()) continue;
			if (c == null) continue;
			if (c.getCigarElements().size() <= 1) continue;
			int hardClipLength = 0;
			List<CigarElement> list = new ArrayList<CigarElement>(c.getCigarElements());
			if (c.getFirstCigarElement().getOperator() == CigarOperator.HARD_CLIP) {
				int length = c.getFirstCigarElement().getLength();
				hardClipLength += length;
				list.set(0, new CigarElement(length, CigarOperator.SOFT_CLIP));
			}
			if (c.getLastCigarElement().getOperator() == CigarOperator.HARD_CLIP) {
				int length = c.getLastCigarElement().getLength();
				hardClipLength += length;
				list.set(list.size() - 1, new CigarElement(length, CigarOperator.SOFT_CLIP));
			}
			if (hardClipLength > 0 && r.getReadLength() + hardClipLength == seq.length) {
				if (r.getReadNegativeStrandFlag()) {
					SequenceUtil.reverseComplement(seq);
					ArrayUtils.reverse(qual);
				}
				r.setReadBases(seq);
				r.setBaseQualities(qual);
				list = CigarUtil.clean(list, false);
				r.setCigar(new Cigar(list));
			}
			if (r.getCigarString().contains("H")) {
				log.warn(String.format("Unable to soften hard clip for %s", r.getReadName()));
			}
		}
	}
	/**
	 * Gets the full read sequence
	 * @param records
	 * @return
	 */
	private static final byte[] getFullSequence(List<SAMRecord> records) {
		if (records == null || records.size() == 0) return null;
		SAMRecord best = records.get(0);
		for (SAMRecord r : records) {
			if (r.getReadBases().length > best.getReadBases().length) {
				best = r;
			}
		}
		byte[] readBases = Arrays.copyOf(best.getReadBases(), best.getReadBases().length);
		if (readBases == null || readBases == SAMRecord.NULL_SEQUENCE) return null;
		if (!best.getReadUnmappedFlag() && best.getReadNegativeStrandFlag()) {
			SequenceUtil.reverseComplement(readBases);
		}
		return readBases;
	}
	/**
	 * Gets the full base quality scores
	 * @param records
	 * @return
	 */
	private static final byte[] getFullBaseQualities(List<SAMRecord> records) {
		if (records == null || records.size() == 0) return null;
		SAMRecord best = records.get(0);
		for (SAMRecord r : records) {
			if (r.getBaseQualities().length > best.getBaseQualities().length) {
				best = r;
			}
		}
		byte[] baseQuals = Arrays.copyOf(best.getBaseQualities(), best.getBaseQualities().length);
		if (baseQuals == null || baseQuals == SAMRecord.NULL_QUALS) return null;
		if (!best.getReadUnmappedFlag() && best.getReadNegativeStrandFlag()) {
			ArrayUtils.reverse(baseQuals);
		}
		return baseQuals;
	}
	/**
	 * The index of segment in the template
	 * @param record SAMRecord
	 * @return The index of segment in the template
	 */
	public static final int getSegmentIndex(SAMRecord record) {
		Integer hiTag = record.getIntegerAttribute(SAMTag.FI.name());
		if (hiTag != null) return hiTag;
		if (record.getReadPairedFlag()) {
			if (record.getFirstOfPairFlag()) return 0;
			if (record.getSecondOfPairFlag()) return 1;
		}
		return 0; // default to the first segment as per sam specs
	}
	/**
	 * Gets the read offset of the first aligned base.
	 * @param r
	 * @return zero based offset from start of read based on fastq sequencing base order
	 */
	public static final int getFirstAlignedBaseReadOffset(SAMRecord r) {
		Cigar c = r.getCigar();
		if (c == null || c.getCigarElements().size() == 0) return -1;
		if (r.getReadNegativeStrandFlag()) {
			return getEndClipLength(c.getCigarElements());
		} else {
			return getStartClipLength(c.getCigarElements());
		}
	}
	/**
	 * Gets the read offset of the final aligned base.
	 * @param r
	 * @return zero based offset of the final aligned read base based on fastq sequencing base order
	 */
	public static final int getLastAlignedBaseReadOffset(SAMRecord r) {
		Cigar c = r.getCigar();
		if (c == null || c.getCigarElements().size() == 0) return -1;
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
	 * Gets a string that is unique for the given alignment of the read.
	 * Note: the identifier is only unique with respect to the alignment of the segment/read
	 * not the entire template/read pair.
	 * @param record record
	 * @return alignment-unique identifier for the given SAMRecord
	 */
	public static String getAlignmentUniqueName(SAMRecord record) {
		return getAlignmentUniqueName(
				record.getReadName(),
				SAMRecordUtil.getSegmentIndex(record),
				record.getReadUnmappedFlag(),
				record.getReferenceName(),
				record.getAlignmentStart(),
				record.getReadNegativeStrandFlag(),
				record.getCigarString());
	}
	public static String getAlignmentUniqueName(String readName, int segmentIndex, boolean unmapped, String referenceName, int alignmentStart, boolean negativeStrand, String cigar) {
		StringBuilder sb = new StringBuilder(readName);
		sb.append(SAMRecordUtil.SUFFIX_SEPERATOR);
		sb.append(segmentIndex);
		if (!unmapped) {
			sb.append(SAMRecordUtil.SUFFIX_SEPERATOR);
			sb.append(referenceName);
			sb.append(SAMRecordUtil.SUFFIX_SEPERATOR);
			sb.append(alignmentStart);
			sb.append(SAMRecordUtil.SUFFIX_SEPERATOR);
			sb.append(negativeStrand ? '-' : '+');
			sb.append(SAMRecordUtil.SUFFIX_SEPERATOR);
			sb.append(cigar);
		}
		;
		return sb.toString();
	}
	/**
	 * Converts a record with mapq below the given mapq threshold to unmapped reads
	 * @param record SAMRecord to convert
	 * @param minMapq minimum MAPQ to avoid unmapping
	 */
	public static SAMRecord lowMapqToUnmapped(SAMRecord record, int minMapq) {
		lowMapqToUnmapped_SA(record, minMapq);
		if (record.getMappingQuality() < minMapq) {
			record.setReadUnmappedFlag(true);
		}
		Integer mateMapq = record.getIntegerAttribute(SAMTag.MQ.name());
		if (mateMapq != null && mateMapq < minMapq) {
			record.setMateUnmappedFlag(true);
		}
		return record;
	}
	private static void lowMapqToUnmapped_SA(SAMRecord record, int minMapq) {
		if (record.getAttribute(SAMTag.SA.name()) == null) return;
		record.setAttribute(SAMTag.SA.name(), ChimericAlignment.getChimericAlignments(record).stream()
				.filter(ca -> ca.mapq >= minMapq)
				.map(ca -> ca.toString())
				.collect(Collectors.joining(";")));
	}
	/**
	 * Converts soft clipped bases that match the reference genome to alignment matches.
	 *  In general, such bases should be included in the alignment but, due to error correction
	 *  of sequencing errors, this does not necessarily hold true for assembly contigs.
	 * @param r record
	 */
	public static void unclipExactReferenceMatches(ReferenceLookup ref, SAMRecord read) {
		if (read.getReadUnmappedFlag()) return;
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
		if (startunclip == 0 && endunclip == 0) return cigar;
		List<CigarElement> list = cigar.getCigarElements();
		list = CigarUtil.clean(list, true);
		if (startunclip > 0) {
			int i = 0;
			while (i < list.size() && list.get(i).getOperator() != CigarOperator.SOFT_CLIP) i++;
			if (i == list.size() || list.get(i).getLength() < startunclip) {
				throw new IllegalArgumentException(String.format("Unable to unclip %d clipped bases from start of %s", startunclip, cigar));
			}
			list.set(i, new CigarElement(list.get(i).getLength() - startunclip, CigarOperator.SOFT_CLIP));
			list.add(i + 1, new CigarElement(startunclip, CigarOperator.MATCH_OR_MISMATCH));
			CigarUtil.clean(list, false);
		}
		if (endunclip > 0) {
			int i = list.size() - 1;
			while (i > 0 && list.get(i).getOperator() != CigarOperator.SOFT_CLIP) i--;
			if (i < 0 || list.get(i).getLength() < endunclip) {
				throw new IllegalArgumentException(String.format("Unable to unclip %d clipped bases from end of %s", endunclip, cigar));
			}
			list.set(i, new CigarElement(list.get(i).getLength() - endunclip, CigarOperator.SOFT_CLIP));
			list.add(i - 1, new CigarElement(endunclip, CigarOperator.MATCH_OR_MISMATCH));
			CigarUtil.clean(list, false);
		}
		return new Cigar(list);
	}
}
