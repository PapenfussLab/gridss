package au.edu.wehi.idsv;

import com.google.common.collect.ImmutableList;




public class RealignedRemoteSAMRecordAssemblyEvidenceTest extends RemoteEvidenceTest {
	@Override
	public RealignedRemoteSAMRecordAssemblyEvidence makeRemote(BreakendSummary bs, String allBases, String realignCigar, boolean realignNegativeStrand) {
		return (RealignedRemoteSAMRecordAssemblyEvidence)makeLocal(bs, allBases, realignCigar, realignNegativeStrand).asRemote();
	}
	@Override
	public RealignedSAMRecordAssemblyEvidence makeLocal(BreakendSummary bs, String allBases, String realignCigar, boolean realignNegativeStrand) {
		RealignedSAMRecordAssemblyEvidence ras = new RealignedSAMRecordAssemblyEvidence(AES(), AssemblyFactory.createAnchoredBreakend(getContext(), AES(), bs.direction, null, bs.referenceIndex, bs.start, anchorLength(allBases, realignCigar), B(allBases), B(allBases)).annotateAssembly().getSAMRecord(),
				ImmutableList.of(realignSAM(bs, allBases, realignCigar, realignNegativeStrand)));
		ras.getSAMRecord().setMappingQuality(21);
		ras.getBackingRecord().setMappingQuality(21);
		return ras;
	}
}
