package au.edu.wehi.idsv;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.NavigableMap;
import java.util.TreeMap;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.tuple.Pair;

import au.edu.wehi.idsv.sam.CigarUtil;
import au.edu.wehi.idsv.sam.SAMRecordUtil;

import com.google.common.collect.ImmutableList;

/**
 * Alignment of a sequence contig that spans multiple segments
 * @author Daniel Cameron
 *
 */
public class CompoundBreakendAlignment {
	private final ProcessingContext processContext;
	private final SAMFileHeader header;
	private final NavigableMap<Integer, SAMRecord> mappedStartOffset = new TreeMap<Integer, SAMRecord>();
	private final BreakendSummary local;
	private final byte[] anchorSequence;
	private final byte[] anchorQual;
	private final byte[] breakendSequence;
	private final byte[] breakendQual;
	public int getBreakpointCount() {
		return mappedStartOffset.size();
	}
	public SAMRecord getSimpleBreakendRealignment() {
		SAMRecord r;
		if (mappedStartOffset.size() == 0) {
			r = getPlaceholderRealignment();
		} else {
			if (local.direction == BreakendDirection.Forward) {
				r = getSoftClipRealignment(mappedStartOffset.firstEntry().getValue());
			} else {
				r = getSoftClipRealignment(mappedStartOffset.lastEntry().getValue());
			}
		}
		return r;
	}
	/**
	 * Gets the set of subsequent breakpoints that the breakend sequence spans
	 * @return
	 */
	public List<Pair<SAMRecord, SAMRecord>> getSubsequentBreakpointAlignmentPairs() {
		if (getBreakpointCount() <= 1) return ImmutableList.of();
		BreakendDirection direction = local.direction;
		byte[] fullSeq;
		byte[] fullQual;
		if (direction == BreakendDirection.Forward) {
			fullSeq = ArrayUtils.addAll(anchorSequence, breakendSequence);
			fullQual = ArrayUtils.addAll(anchorQual, breakendQual);
		} else {
			fullSeq = ArrayUtils.addAll(breakendSequence, anchorSequence);
			fullQual = ArrayUtils.addAll(breakendQual, anchorQual);
		}
		List<Pair<SAMRecord, SAMRecord>> output = new ArrayList<Pair<SAMRecord,SAMRecord>>(getBreakpointCount() - 1);
		SAMRecord lastEntry = null;
		int offsetDueToAnchor = (local.direction == BreakendDirection.Forward ? anchorSequence.length : 0);
		for (SAMRecord entry : mappedStartOffset.values()) {
			if (lastEntry != null) {
				SAMRecord left = asFullSequence(fullSeq, fullQual, lastEntry, getReadOffset(lastEntry) + offsetDueToAnchor);
				SAMRecord right = asFullSequence(fullSeq, fullQual, entry, getReadOffset(entry) + offsetDueToAnchor);
				if (local.direction == BreakendDirection.Forward) {
					int anchorBases = fullSeq.length - (left.getReadNegativeStrandFlag() ? SAMRecordUtil.getStartSoftClipLength(left) : SAMRecordUtil.getEndSoftClipLength(left));
					trimSoftClip(right.getReadNegativeStrandFlag() ? 0 : anchorBases, right.getReadNegativeStrandFlag() ? anchorBases : 0, right);
				} else {
					int anchorBases = fullSeq.length - (right.getReadNegativeStrandFlag() ? SAMRecordUtil.getEndSoftClipLength(right) : SAMRecordUtil.getStartSoftClipLength(right));
					trimSoftClip(left.getReadNegativeStrandFlag() ? anchorBases : 0, left.getReadNegativeStrandFlag() ? 0 : anchorBases, left);
				}
				output.add(Pair.of(left, right));
			}
			lastEntry = entry;
		}
		return output;
	}
	private void trimSoftClip(int start, int end, SAMRecord read) {
		int length = read.getReadLength();
		read.setReadBases(Arrays.copyOfRange(read.getReadBases(), start, length - end));
		read.setBaseQualities(Arrays.copyOfRange(read.getBaseQualities(), start, length - end));
		List<CigarElement> list = new ArrayList<CigarElement>(read.getCigar().getCigarElements());
		if (start > 0) {
			CigarElement e = list.get(0);
			assert(e.getOperator() == CigarOperator.SOFT_CLIP);
			assert(e.getLength() >= start);
			list.set(0, new CigarElement(e.getLength() - start, CigarOperator.SOFT_CLIP));
		}
		if (end > 0) {
			CigarElement e = list.get(list.size() - 1);
			assert(e.getOperator() == CigarOperator.SOFT_CLIP);
			assert(e.getLength() >= start);
			list.set(list.size() - 1, new CigarElement(e.getLength() - end, CigarOperator.SOFT_CLIP));
		}
		if (list.get(0).getOperator() == CigarOperator.SOFT_CLIP && list.get(0).getLength() == 0) {
			list.remove(0);
		}
		if (list.get(list.size() - 1).getOperator() == CigarOperator.SOFT_CLIP && list.get(list.size() - 1).getLength() == 0) {
			list.remove(list.size() - 1);
		}
		read.setCigar(new Cigar(list));
	}
	private SAMRecord asFullSequence(byte[] fullSeq, byte[] fullQual, SAMRecord realign, int readStartOffset) {
		SAMRecord r = SAMRecordUtil.clone(realign);
		assert(r.getReadLength() == r.getCigar().getReadLength()); // sanity check cigar matches read base count
		int length = r.getReadLength() - SAMRecordUtil.getStartSoftClipLength(r) - SAMRecordUtil.getEndSoftClipLength(r);
		int startPad;
		int endPad;
		byte[] seq = fullSeq;
		byte[] qual = fullQual;
		if (r.getReadNegativeStrandFlag()) {
			seq = Arrays.copyOf(seq, seq.length);
			qual = Arrays.copyOf(qual, qual.length);
			SequenceUtil.reverseComplement(seq);
			SequenceUtil.reverseComplement(qual);
			endPad = readStartOffset + SAMRecordUtil.getEndSoftClipLength(r);
			startPad = fullSeq.length - length - endPad;
		} else {
			startPad = readStartOffset + SAMRecordUtil.getStartSoftClipLength(r);
			endPad = fullSeq.length - length - startPad;
		}
		List<CigarElement> anchorCigarList = new ArrayList<CigarElement>(r.getCigar().getCigarElements());
		CigarUtil.trimClipping(anchorCigarList);
		CigarUtil.addStartSoftClip(anchorCigarList, startPad);
		CigarUtil.addEndSoftClip(anchorCigarList, endPad);
		r.setReadBases(seq);
		r.setBaseQualities(qual);
		r.setCigar(new Cigar(anchorCigarList));
		assert(r.getCigar().getReadLength() == fullSeq.length);
		return r;
	}
	private SAMRecord getSoftClipRealignment(SAMRecord realign) {
		if (realign.getReadLength() == breakendSequence.length) return realign;
		return asFullSequence(breakendSequence, breakendQual, realign, getReadOffset(realign));
	}
	private SAMRecord getPlaceholderRealignment() {
		SAMRecord placeholder = new SAMRecord(header);
		placeholder.setReadUnmappedFlag(true);
		placeholder.setReadBases(breakendSequence);
		placeholder.setBaseQualities(breakendQual);
		placeholder.setReadNegativeStrandFlag(false);
		placeholder.setMappingQuality(SAMRecord.NO_MAPPING_QUALITY);
		return placeholder;
	}
	private boolean isFiltered(SAMRecord realignment) {
		return realignment == null ||
				realignment.getReadUnmappedFlag() ||
				realignment.getMappingQuality() < processContext.getRealignmentParameters().mapqUniqueThreshold;
	}
	public CompoundBreakendAlignment(
			ProcessingContext processContext,
			SAMFileHeader header,
			BreakendSummary local,
			byte[] anchorSequence,
			byte[] anchorQual,
			byte[] breakendSequence,
			byte[] breakendQual,
			Collection<SAMRecord> realignments) {
		this.processContext = processContext;
		this.header = header;
		this.local = local;
		this.anchorSequence = anchorSequence;
		this.anchorQual = anchorQual;
		this.breakendSequence = breakendSequence;
		this.breakendQual = breakendQual;
		if (realignments != null) {
			for (SAMRecord r : realignments) {
				if (r != null && !isFiltered(r)) {
					mappedStartOffset.put(getMappedOffset(r), r);
				}
			}
		}
	}
	/**
	 * Gets the breakend offset of the first base
	 * @param r realignment
	 * @return breakend offset of the first read base
	 */
	private int getReadOffset(SAMRecord r) {
		if (r.getReadLength() == breakendSequence.length) return 0;
		return BreakpointFastqEncoding.getEncodedBreakendOffset(r.getReadName());
	}
	/**
	 * Gets the breakend offset of the first mapped base
	 * @param r realignment
	 * @return breakend offset of the first mapped read base
	 */
	private int getMappedOffset(SAMRecord r) {
		int offset = getReadOffset(r);
		if (r.getReadNegativeStrandFlag()) {
			offset += SAMRecordUtil.getEndSoftClipLength(r);
		} else {
			offset += SAMRecordUtil.getStartSoftClipLength(r);
		}
		return offset;
	}
}
