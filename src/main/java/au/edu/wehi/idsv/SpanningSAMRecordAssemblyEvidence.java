package au.edu.wehi.idsv;

import java.util.List;

import htsjdk.samtools.SAMRecord;
import au.edu.wehi.idsv.sam.SplitIndel;

import com.google.common.collect.ImmutableList;

/**
 * Assembly spanning a small indel
 */
public class SpanningSAMRecordAssemblyEvidence extends RealignedSAMRecordAssemblyEvidence {
	private final SAMRecordAssemblyEvidence parent;
	private final String evidenceID;
	private SpanningSAMRecordAssemblyEvidence remote;
	/**
	 * Creates the forward breakpoint of the event spanned by the assembly
	 * @param source
	 * @param assembly
	 */
	public SpanningSAMRecordAssemblyEvidence(SAMRecordAssemblyEvidence parent, SplitIndel indel, int indelOffset) {
		this(parent, indel.leftAnchored, indel.leftRealigned, indelOffset);
		this.remote = new SpanningSAMRecordAssemblyEvidence(parent, indel.rightAnchored, indel.rightRealigned, indelOffset);
		this.remote.remote = this;
	}
	private SpanningSAMRecordAssemblyEvidence(SAMRecordAssemblyEvidence parent, SAMRecord anchor, SAMRecord realign, int indelOffset) {
		super(parent.getEvidenceSource(), anchor, ImmutableList.of(realign));
		this.parent = parent;
		this.evidenceID = String.format("%s_%s%d", parent.getEvidenceID(), getBreakendSummary().direction.toChar(), indelOffset);
		// sanity check mapping quality
		assert(anchor.getMappingQuality() >= parent.getEvidenceSource().getContext().getConfig().minMapq);
	}
	public SAMRecordAssemblyEvidence getParentAssembly() {
		return parent;
	}
	@Override
	public SpanningSAMRecordAssemblyEvidence asRemote() {
		return remote;
	}
	@Override
	public String getEvidenceID() {
		return evidenceID;
	}
	@Override
	public SAMRecord getBackingRecord() {
		return parent.getSAMRecord();
	}
	@Override
	public SpanningSAMRecordAssemblyEvidence annotateAssembly() {
		throw new UnsupportedOperationException();
	}
	@Override
	public List<SpanningSAMRecordAssemblyEvidence> getSpannedIndels() {
		throw new UnsupportedOperationException();
	}
	@Override
	public float getBreakendQual() {
		return super.getBreakendQual() / 2;
	}
	@Override
	public float getBreakpointQual() {
		return super.getBreakpointQual() / 2;
	}
}
