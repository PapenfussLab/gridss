package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;

import java.util.Collection;

import au.edu.wehi.idsv.picard.ReferenceLookup;
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
	private final RangeMap<Integer, SAMRecord> map = TreeRangeMap.create();
	private final BreakendSummary local;
	private final byte[] anchoredSequence;
	public int getBreakpointCount() {
		return map.asMapOfRanges().size();
	}
	public BreakendSummary getPrimaryBreakend() {
	}
	public SAMRecord getSimpleBreakendRealignment(SAMRecord anchor) {
		SAMRecord r;
		if (map.asMapOfRanges().size() == 0) {
			r = getPlaceholderRealignment(anchor, anchoredSequence, );
		} else if (map.asMapOfRanges().size() == 1) {
			// technically not correct as the first alignment could be multimapping
			r = map.asMapOfRanges().values().iterator().next();
		}
		int breakendLength = map.span().upperEndpoint() - map.span().lowerEndpoint();
		
		SAMRecordUtil.pairReads(anchor, r);
		return r;
	}
	private static SAMRecord getPlaceholderRealignment(SAMRecord anchor, byte[] breakendSequence, byte[] breakendQual) {
		SAMRecord placeholder = new SAMRecord(anchor.getHeader());
		placeholder.setReadUnmappedFlag(true);
		placeholder.setReadBases(breakendSequence);
		placeholder.setBaseQualities(breakendQual);
		placeholder.setReadNegativeStrandFlag(false);
		placeholder.setMappingQuality(0);
		return placeholder;
	}
	private static void pairReads(SAMRecord anchor, SAMRecord realignment) {
		if (realignment.getReadUnmappedFlag()) {
			// SAMv1 S2.4
			realignment.setReferenceIndex(anchor.getReferenceIndex());
			realignment.setAlignmentStart(anchor.getAlignmentStart());
		}
		SAMRecordUtil.pairReads(anchor, realignment);
	}
	public CompoundBreakendAlignment(
			RealignmentParameters rp,
			BreakendSummary local,
			byte[] anchoredSequence,
			Collection<SAMRecord> realignments) {
		this.local = local;
		this.anchoredSequence = anchoredSequence;
		for (SAMRecord r : realignments) {
			if (r != null && !r.getReadUnmappedFlag() && r.getMappingQuality() >= rp.mapqUniqueThreshold) {
				int offset = BreakpointFastqEncoding.getEncodedReadOffset(r.getReadName());
				map.put(Range.closedOpen(offset, offset + r.getReadLength()), r);
			} else {
				// ignore alignments that aren't unique enough
			}
		}
	}
}
