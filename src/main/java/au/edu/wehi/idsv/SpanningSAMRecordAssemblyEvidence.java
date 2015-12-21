package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;

import com.google.common.collect.ImmutableList;

/**
 * Assembly spanning a small indel
 */
public class SpanningSAMRecordAssemblyEvidence extends RealignedSAMRecordAssemblyEvidence {
	private final SAMRecordAssemblyEvidence parent;
	private final String evidenceID;
	/**
	 * Creates the forward breakpoint of the event spanned by the assembly
	 * @param source
	 * @param assembly
	 */
	public SpanningSAMRecordAssemblyEvidence(SAMRecordAssemblyEvidence parent, SAMRecord splitAnchor, SAMRecord splitRealign, int indelOffset) {
		super(parent.getEvidenceSource(), splitAnchor, ImmutableList.of(splitRealign));
		this.parent = parent;
		this.evidenceID = getSAMRecord().getReadName() + getBreakendSummary().direction.toChar();
	}
	public SAMRecordAssemblyEvidence getParentAssembly() {
		return parent;
	}
	@Override
	public SAMRecord getBackingRecord() {
		throw new UnsupportedOperationException();
	}
	@Override
	public SpanningSAMRecordAssemblyEvidence asRemote() {
		throw new UnsupportedOperationException();
	}
	@Override
	public String getEvidenceID() {
		assert(parent != null);
		return evidenceID;
	}
	@Override
	public boolean isSpanningAssembly() {
		return true;
	}
	@Override
	public SpanningSAMRecordAssemblyEvidence annotateAssembly() {
		throw new UnsupportedOperationException();
	}
}
