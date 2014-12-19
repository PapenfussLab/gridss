package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordComparator;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.sam.SamTags;

import com.google.common.collect.ComparisonChain;

public class RealignedSoftClipEvidence extends SoftClipEvidence implements DirectedBreakpoint {
	private final RealignedBreakpoint rbp;
	private final SAMRecord realigned;
	protected RealignedSoftClipEvidence(ProcessingContext processContext, SAMEvidenceSource source, BreakendDirection direction, SAMRecord record, SAMRecord realigned) {
		super(processContext, source, direction, SAMRecordUtil.clone(record)); // need a copy since we're changing attributes to different values for the forward and backward records 
		this.realigned = realigned;
		int pos = direction == BreakendDirection.Forward ? record.getAlignmentEnd() : record.getAlignmentStart();
		BreakendSummary local = new BreakendSummary(record.getReferenceIndex(), direction, pos, pos);
		this.rbp = new RealignedBreakpoint(processContext, local, record.getReadBases(), realigned);
		setPositionAttributes();
	}
	private void setPositionAttributes() {
		SAMRecord sc = getSAMRecord();
		sc.setAttribute(SamTags.REALIGNMENT_REFERENCE_INDEX, realigned.getReferenceIndex());
		sc.setAttribute(SamTags.REALIGNMENT_POSITION, realigned.getAlignmentStart());
		/*
		realigned.setAttribute(SamTags.REALIGNMENT_REFERENCE_INDEX, sc.getReferenceIndex());
		realigned.setAttribute(SamTags.REALIGNMENT_POSITION, sc.getAlignmentStart());
		realigned.setReadPairedFlag(true);
		realigned.setMateUnmappedFlag(false);
		realigned.setMateReferenceIndex(sc.getReferenceIndex());
		realigned.setMateAlignmentStart(sc.getAlignmentStart());
		*/
	}
	public SAMRecord getRealignedSAMRecord() {
		return realigned;
	}
	@Override
	public BreakpointSummary getBreakendSummary() {
		return rbp.getBreakpointSummary();
	}
	public String getUntemplatedSequence() {
		return rbp.getInsertedSequence();
	}
	@Override
	public int getRemoteMapq() {
		return realigned.getMappingQuality();
	}
	@Override
	public int getRemoteBaseLength() {
		return realigned.getReadLength();
	}
	@Override
	public int getRemoteBaseCount() {
		return getRemoteBaseLength();
	}
	@Override
	public int getRemoteMaxBaseQual() {
		return SAMRecordUtil.getMaxReferenceBaseQual(realigned);
	}
	@Override
	public int getRemoteTotalBaseQual() {
		return SAMRecordUtil.getTotalReferenceBaseQual(realigned);
	}
	/**
	 * Comparator for sorting SAMRecords by coordinate of the matched realignment.
	 *
	 */
	public static class RealignmentCoordinateComparator implements SAMRecordComparator {
		public int compare(final SAMRecord arg0, final SAMRecord arg1) {
	    	return ComparisonChain.start()
			        .compare(arg0.getIntegerAttribute(SamTags.REALIGNMENT_REFERENCE_INDEX), arg1.getIntegerAttribute(SamTags.REALIGNMENT_REFERENCE_INDEX))
			        .compare(arg0.getIntegerAttribute(SamTags.REALIGNMENT_POSITION), arg1.getIntegerAttribute(SamTags.REALIGNMENT_POSITION))
			        .compare(arg0.getReadName(), arg1.getReadName())
			        .compare(arg0.getFlags(), arg1.getFlags())
			        .result();
	    }
	    /**
	     * @return negative if samRecord1 < samRecord2,  0 if equal, else positive
	     */
	    @Override
		public int fileOrderCompare(final SAMRecord arg0, final SAMRecord arg1) {
	    	return ComparisonChain.start()
		        .compare(arg0.getIntegerAttribute(SamTags.REALIGNMENT_REFERENCE_INDEX), arg1.getIntegerAttribute(SamTags.REALIGNMENT_REFERENCE_INDEX))
		        .compare(arg0.getIntegerAttribute(SamTags.REALIGNMENT_POSITION), arg1.getIntegerAttribute(SamTags.REALIGNMENT_POSITION))
		        .result();
	    }
	}
	@Override
	public double getBreakpointQual() {
		return scPhred(getEvidenceSource(), getSoftClipLength(), getLocalMapq(), getRemoteMapq());
	}
}
