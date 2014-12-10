package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

import org.junit.Test;


public abstract class RemoteEvidenceTest extends TestHelper  {
	protected int anchorLength(final String allBases, final String realignCigar) {
		return allBases.length() - (new SAMRecord(null) {{ setCigarString(realignCigar); }}).getCigar().getReadLength();
	}
	protected SAMRecord realignSAM(final BreakendSummary bs, final String allBases, final String realignCigar, final boolean realignNegativeStrand) {
		return new SAMRecord(getContext().getBasicSamHeader()) {{
			setReferenceIndex(1);
			setAlignmentStart(100);
			setCigarString(realignCigar);
			setReadNegativeStrandFlag(realignNegativeStrand);
			String realignBases;
			if (bs.direction == FWD) {
				realignBases = allBases.substring(anchorLength(allBases, realignCigar));
			} else {
				realignBases = allBases.substring(0, allBases.length() - anchorLength(allBases, realignCigar));
			}
			if (realignNegativeStrand) {
				realignBases = SequenceUtil.reverseComplement(realignBases);
			}
			setReadBases(B(realignBases));
			setMappingQuality(21);
			}};
	}
	public abstract RemoteEvidence makeRemote(BreakendSummary bs, final String allBases, final String realignCigar, final boolean realignNegativeStrand);
	public abstract DirectedBreakpoint makeLocal(BreakendSummary bs, final String allBases, final String realignCigar, final boolean realignNegativeStrand);
	private RemoteEvidence R() { return makeRemote(new BreakendSummary(0, FWD, 1, 1), "TTTTTTTTTT", "4M", true); } 
	private DirectedBreakpoint L() { return makeLocal(new BreakendSummary(0, FWD, 1, 1), "TTTTTTTTTT", "4M", true); }
	@Test
	public void getUntemplatedSequence_should_return_untemplated_sequence_of_local_side_of_breakend() throws CloneNotSupportedException {
		// which is just whatever the soft clip sequences of the aligned read is (as the mapper handles the reverse complimenting)
		assertEquals("GTN", makeRemote(new BreakendSummary(0, FWD, 1, 1), "NGTNCA", "3S2M", false).getUntemplatedSequence());
		assertEquals(SequenceUtil.reverseComplement("GTN"), makeRemote(new BreakendSummary(0, FWD, 1, 1), "NGTNCA", "2M3S", true).getUntemplatedSequence());
		assertEquals("NCA", makeRemote(new BreakendSummary(0, BWD, 1, 1), "GTNCAN", "2M3S", false).getUntemplatedSequence());
		assertEquals(SequenceUtil.reverseComplement("NCA"), makeRemote(new BreakendSummary(0, BWD, 1, 1), "GTNCAN", "3S2M", true).getUntemplatedSequence());
	}
	@Test
	public void getBreakendSequence_should_return_sequence_of_local_side() throws CloneNotSupportedException {
		// TCGTNCA> <GTNCA
		assertEquals("TCGTN", S(makeRemote(new BreakendSummary(0, FWD, 1, 1), "TCGTNCA", "3S2M", false).getBreakendSequence()));
		// TCGTNCA>   RC>
		assertEquals(SequenceUtil.reverseComplement("TCGTN"), S(makeRemote(new BreakendSummary(0, FWD, 1, 1), "TCGTNCA", "2M3S", true).getBreakendSequence()));
		// TCGTNCA>  <TCGTNCA
		assertEquals("GTNCA", S(makeRemote(new BreakendSummary(0, BWD, 1, 1), "TCGTNCA", "2M3S", false).getBreakendSequence()));
		
		assertEquals(SequenceUtil.reverseComplement("GTNCA"), S(makeRemote(new BreakendSummary(0, BWD, 1, 1), "TCGTNCA", "3S2M", true).getBreakendSequence()));
	}
	@Test
	public void should_swap_breakpoint() {
		assertEquals(L().getBreakendSummary().remoteBreakpoint(), R().getBreakendSummary());
	}
	@Test
	public void should_swap_remote_local_fields() {
		assertEquals(L().getLocalMapq(), R().getRemoteMapq());
		assertEquals(L().getLocalBaseLength(), R().getRemoteBaseLength());
		assertEquals(L().getLocalMaxBaseQual(), R().getRemoteMaxBaseQual());
		assertEquals(L().getLocalTotalBaseQual(), R().getRemoteTotalBaseQual());
		
		assertEquals(R().getLocalMapq(), L().getRemoteMapq());
		assertEquals(R().getLocalBaseLength(), L().getRemoteBaseLength());
		assertEquals(R().getLocalMaxBaseQual(), L().getRemoteMaxBaseQual());
		assertEquals(R().getLocalTotalBaseQual(), L().getRemoteTotalBaseQual());
	}
}
