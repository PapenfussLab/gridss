package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;


public class RealignedRemoteSoftClipEvidenceTest extends RemoteEvidenceTest {
	@Override
	public RemoteEvidence makeRemote(BreakendSummary bs, String allBases, String realignCigar, boolean realignNegativeStrand) {
		try {
			return new RealignedRemoteSoftClipEvidence(getContext(), SES(), bs.direction, makeSC(bs, allBases, realignCigar), realignSAM(bs, allBases, realignCigar, realignNegativeStrand));
		} catch (CloneNotSupportedException e) {
			throw new RuntimeException();
		}
	}
	@Override
	public DirectedBreakpoint makeLocal(BreakendSummary bs, String allBases, String realignCigar, boolean realignNegativeStrand) {
		try {
			return new RealignedSoftClipEvidence(getContext(), SES(), bs.direction, makeSC(bs, allBases, realignCigar), realignSAM(bs, allBases, realignCigar, realignNegativeStrand));
		} catch (CloneNotSupportedException e) {
			throw new RuntimeException();
		}
	}
	private SAMRecord makeSC(final BreakendSummary bs, final String allBases, final String realignCigar) {
		return new SAMRecord(getContext().getBasicSamHeader()) {{
			setReferenceIndex(bs.referenceIndex);
			setReadBases(B(allBases));
			setMappingQuality(40);
			int anchorLen = anchorLength(allBases, realignCigar);
			if (bs.direction == FWD) {
				setCigarString(String.format("%dM%dS", anchorLen, allBases.length() - anchorLen));
				setAlignmentStart(bs.start - anchorLen + 1);
			} else {
				setCigarString(String.format("%dS%dM", allBases.length() - anchorLen, anchorLen));
				setAlignmentStart(bs.start);
			}
		}};
	}
}
