package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;


public class RealignedRemoteSoftClipEvidenceTest extends RemoteEvidenceTest {
	@Override
	public RealignedRemoteSoftClipEvidence makeRemote(BreakendSummary bs, String allBases, String realignCigar, boolean realignNegativeStrand) {
		return makeLocal(bs, allBases, realignCigar, realignNegativeStrand).asRemote();
	}
	@Override
	public RealignedSoftClipEvidence makeLocal(BreakendSummary bs, String allBases, String realignCigar, boolean realignNegativeStrand) {
		return new RealignedSoftClipEvidence(SES(), bs.direction, makeSC(bs, allBases, realignCigar), realignSAM(bs, allBases, realignCigar, realignNegativeStrand));
	}
	private SAMRecord makeSC(final BreakendSummary bs, final String allBases, final String realignCigar) {
		return new SAMRecord(getContext().getBasicSamHeader()) {{
			setReadName("readName");
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
	@Test
	public void evidenceID_should_correspond_to_underlying_evidenceID() {
		SAMRecord read = Read(0, 1, "1M1S");
		read.setReadPairedFlag(true);
		read.setFirstOfPairFlag(true);
		read.setMateUnmappedFlag(true);
		RealignedSoftClipEvidence realigned = (RealignedSoftClipEvidence) SoftClipEvidence.create(SES(), FWD, read, Read(0, 5, "1M"));
		assertEquals("R" + realigned.getEvidenceID(), realigned.asRemote().getEvidenceID());
		
		read.setFirstOfPairFlag(false);
		read.setSecondOfPairFlag(true);
		realigned = (RealignedSoftClipEvidence) SoftClipEvidence.create(SES(), FWD, read, Read(0, 5, "1M"));
		assertEquals("R" + realigned.getEvidenceID(), realigned.asRemote().getEvidenceID());
	}
}
