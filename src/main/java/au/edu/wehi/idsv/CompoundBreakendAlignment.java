package au.edu.wehi.idsv;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Map.Entry;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.tuple.Pair;

import au.edu.wehi.idsv.sam.CigarUtil;
import au.edu.wehi.idsv.sam.SAMRecordUtil;

import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.TreeRangeMap;

/**
 * Alignment of a sequence contig that spans multiple segments
 * @author cameron.d
 *
 */
public class CompoundBreakendAlignment {
	private final ProcessingContext processContext;
	private final SAMFileHeader header;
	private final RangeMap<Integer, SAMRecord> map = TreeRangeMap.create();
	private final BreakendSummary local;
	private final byte[] anchorSequence;
	private final byte[] anchorQual;
	private final byte[] breakendSequence;
	private final byte[] breakendQual;
	public int getBreakpointCount() {
		return map.asMapOfRanges().size();
	}
	public SAMRecord getSimpleBreakendRealignment() {
		SAMRecord r;
		if (map.asMapOfRanges().size() == 0) {
			r = getPlaceholderRealignment();
		} else {
			if (local.direction == BreakendDirection.Forward) {
				r = getSoftClipRealignment(map.getEntry(map.span().lowerEndpoint()));	
			} else {
				r = getSoftClipRealignment(map.getEntry(map.span().upperEndpoint() - 1));
			}
		}
		return r;
	}
	/**
	 * Gets the set of subsequent breakpoints that the breakend sequence spans
	 * @return
	 */
	public List<Pair<SAMRecord, SAMRecord>> getSubsequentBreakpointAlignmentPairs() {
		BreakendDirection direction = local.direction;
		int offsetShift = 0;
		byte[] fullSeq;
		byte[] fullQual;
		if (direction == BreakendDirection.Forward) {
			fullSeq = ArrayUtils.addAll(anchorSequence, breakendSequence);
			fullQual = ArrayUtils.addAll(anchorQual, breakendQual);
			offsetShift = anchorSequence.length;
		} else {
			fullSeq = ArrayUtils.addAll(breakendSequence, anchorSequence);
			fullQual = ArrayUtils.addAll(breakendQual, anchorQual);
		}
		List<Pair<SAMRecord, SAMRecord>> output = new ArrayList<Pair<SAMRecord,SAMRecord>>(getBreakpointCount() - 1);
		Entry<Range<Integer>, SAMRecord> lastEntry = null;
		for (Entry<Range<Integer>, SAMRecord> entry : map.asMapOfRanges().entrySet()) {
			if (lastEntry != null) {
				// process entries
				// TODO: process these
				throw new RuntimeException("NYI");
			}
			lastEntry = entry;
		}
		return output;
	}
	private SAMRecord getSoftClipRealignment(Entry<Range<Integer>, SAMRecord> entry) {
		return getSoftClipRealignment(entry.getKey().lowerEndpoint(), entry.getValue());
	}
	/**
	 * Converts the given partial realignment into an equivalent realignment
	 * record for the full breakend sequence 
	 * @param offsetStart
	 * @param offsetEnd
	 * @param r realignment record containing a subset of the breakend sequence
	 * @return
	 */
	private SAMRecord getSoftClipRealignment(int offsetStart, SAMRecord r) {
		if (r.getReadLength() == breakendSequence.length) return r;
		r = SAMRecordUtil.clone(r);
		assert(!r.getReadUnmappedFlag());
		List<CigarElement> cigar = new ArrayList<CigarElement>(r.getCigar().getCigarElements());
		int expandStartSoftClipBy = offsetStart;
		int expandEndSoftClipBy = breakendSequence.length - r.getReadLength();
		if (r.getReadNegativeStrandFlag()) {
			CigarUtil.addStartSoftClip(cigar, expandEndSoftClipBy);
			CigarUtil.addEndSoftClip(cigar, expandStartSoftClipBy);
		} else {
			CigarUtil.addStartSoftClip(cigar, expandStartSoftClipBy);
			CigarUtil.addEndSoftClip(cigar, expandEndSoftClipBy);
		}
		byte[] seq = breakendSequence;
		byte[] qual = breakendQual;
		if (r.getReadNegativeStrandFlag()) {
			seq = Arrays.copyOf(seq, seq.length);
			SequenceUtil.reverseComplement(seq);
			qual = Arrays.copyOf(qual, qual.length);
			ArrayUtils.reverse(qual);
		}
		r.setReadBases(seq);
		r.setBaseQualities(qual);
		r.setCigar(new Cigar(cigar));
		return r;
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
			if (realignments.size() == 1) {
				// simple realignment
				SAMRecord r = realignments.iterator().next();
				if (r != null && !isFiltered(r)) {
					map.put(Range.closedOpen(0, r.getReadLength()), r);
				}
			} else {
				// compound realignment
				for (SAMRecord r : realignments) {
					if (r != null && !isFiltered(r)) {
						int offset = BreakpointFastqEncoding.getEncodedReadOffset(r.getReadName());
						map.put(Range.closedOpen(offset, offset + r.getReadLength()), r);
					}
				}
			}
		}
	}
}
