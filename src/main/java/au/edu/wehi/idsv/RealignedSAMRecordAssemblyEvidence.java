package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;

import java.util.List;

import au.edu.wehi.idsv.sam.SAMRecordUtil;

public class RealignedSAMRecordAssemblyEvidence extends SAMRecordAssemblyEvidence implements DirectedBreakpoint {
	private final RealignedBreakpoint rbp;
	public RealignedSAMRecordAssemblyEvidence(
			AssemblyEvidenceSource source,
			SAMRecordAssemblyEvidence assembly,
			List<SAMRecord> realigned) {
		super(source, assembly, realigned);
		this.rbp = RealignedBreakpoint.create(source.getContext().getReference(), source.getContext().getDictionary(), super.getBreakendSummary(), super.getAssemblyAnchorSequence(), getRemoteSAMRecord());
	}
	public RealignedSAMRecordAssemblyEvidence(
			AssemblyEvidenceSource source,
			SAMRecord assembly,
			List<SAMRecord> realigned) {
		super(source, assembly, realigned);
		this.rbp = RealignedBreakpoint.create(source.getContext().getReference(), source.getContext().getDictionary(), super.getBreakendSummary(), super.getAssemblyAnchorSequence(), getRemoteSAMRecord());
	}
	@Override
	public boolean isReferenceAssembly() {
		return super.isReferenceAssembly()
			|| (getBreakendSummary().couldBeReferenceAllele() && getUntemplatedSequence().length() == 0);
	}
	@Override
	public BreakpointSummary getBreakendSummary() {
		return rbp.getBreakpointSummary();
	}
	@Override
	public String getHomologySequence() {
		return rbp.getHomologySequence();
	}
	@Override
	public int getHomologyAnchoredBaseCount() {
		return rbp.getHomologyBaseCountIncludedLocally();
	}
	@Override
	public String getUntemplatedSequence() {
		return rbp.getInsertedSequence();
	}
	@Override
	public int getRemoteMapq() {
		return getRemoteSAMRecord().getMappingQuality();
	}
	@Override
	public int getRemoteBaseLength() {
		return getRemoteSAMRecord().getReadLength();
	}
	@Override
	public int getRemoteBaseCount() {
		return getRemoteBaseLength();
	}
	@Override
	public int getRemoteMaxBaseQual() {
		return SAMRecordUtil.getMaxReferenceBaseQual(getRemoteSAMRecord());
	}
	@Override
	public int getRemoteTotalBaseQual() {
		return SAMRecordUtil.getTotalReferenceBaseQual(getRemoteSAMRecord());
	}
	public float getBreakpointQual() {
		int rp = getAssemblySupportCountReadPair();
		double rpq = getAssemblySupportReadPairQualityScore();
		int sc = getAssemblySupportCountSoftClip();
		double scq =  getAssemblySupportSoftClipQualityScore();
		if (getEvidenceSource().getContext().getAssemblyParameters().excludeNonSupportingEvidence) {
			rp -= getAssemblyNonSupportingReadPairCount();
			rpq -= getAssemblyNonSupportingReadPairQualityScore();
			sc -= getAssemblyNonSupportingSoftClipCount();
			scq -= getAssemblyNonSupportingSoftClipQualityScore();
		}
		return (float)getEvidenceSource().getContext().getConfig().getScoring().getModel().scoreAssembly(
				rp, rpq,
				sc, scq,
				getLocalMapq(),
				getRemoteMapq());
	}
	/**
	 * Swaps the view of the evidence to the remote breakend 
	 * @return
	 */
	public RealignedSAMRecordAssemblyEvidence asRemote() {
		return new RealignedRemoteSAMRecordAssemblyEvidence(this);
	}
}
