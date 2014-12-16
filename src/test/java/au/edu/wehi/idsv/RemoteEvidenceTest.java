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
	private RemoteEvidence R() { return makeRemote(new BreakendSummary(0, FWD, 100, 100), "GGGGGGGGGG", "4M", true); } 
	private DirectedBreakpoint L() { return makeLocal(new BreakendSummary(0, FWD, 100, 100), "GGGGGGGGGG", "4M", true); }
	@Test
	public void makeRemote_sanityCheck() {
		assertEquals(new BreakendSummary(0, FWD, 100, 100), R().getBreakendSummary().remoteBreakend());
	}
	@Test
	public void makeLocal_sanityCheck() {
		assertEquals(new BreakendSummary(0, FWD, 100, 100), L().getBreakendSummary().localBreakend());
	}
	@Test
	public void getUntemplatedSequence_should_return_untemplated_sequence_of_remote_side_of_breakend() throws CloneNotSupportedException {
		// which is just whatever the soft clip sequences of the aligned read is (as the mapper handles the reverse complimenting)
		assertEquals("GTN", makeRemote(new BreakendSummary(0, FWD, 100, 100), "NGTNCA", "3S2M", false).getUntemplatedSequence());
		assertEquals(SequenceUtil.reverseComplement("GTN"), makeRemote(new BreakendSummary(0, FWD, 100, 100), "NGTNCA", "2M3S", true).getUntemplatedSequence());
		assertEquals("NCA", makeRemote(new BreakendSummary(0, BWD, 100, 100), "GTNCAN", "2M3S", false).getUntemplatedSequence());
		assertEquals(SequenceUtil.reverseComplement("NCA"), makeRemote(new BreakendSummary(0, BWD, 100, 100), "GTNCAN", "3S2M", true).getUntemplatedSequence());
	}
	@Test
	public void getBreakendSequence_should_return_sequence_from_remote_perspective() throws CloneNotSupportedException {
		// TCGTNCA> <GTNCA
		assertEquals("TCGTN", S(makeRemote(new BreakendSummary(0, FWD, 100, 100), "TCGTNCA", "3S2M", false).getBreakendSequence()));
		// TCGTNCA>   RC>
		assertEquals(SequenceUtil.reverseComplement("TCGTN"), S(makeRemote(new BreakendSummary(0, FWD, 100, 100), "TCGTNCA", "2M3S", true).getBreakendSequence()));
		// TCGTNCA>  <TCGTNCA
		assertEquals("GTNCA", S(makeRemote(new BreakendSummary(0, BWD, 100, 100), "TCGTNCA", "2M3S", false).getBreakendSequence()));
		
		assertEquals(SequenceUtil.reverseComplement("GTNCA"), S(makeRemote(new BreakendSummary(0, BWD, 100, 100), "TCGTNCA", "3S2M", true).getBreakendSequence()));
	}
	@Test
	public void should_swap_breakpoint() {
		assertEquals(L().getBreakendSummary().remoteBreakpoint(), R().getBreakendSummary());
	}
	@Test
	public void should_swap_remote_local_fields() {
		DirectedBreakpoint l = L();
		DirectedBreakpoint r = R();
		assertEquals(l.getLocalMapq(), r.getRemoteMapq());
		assertEquals(l.getLocalBaseLength(), r.getRemoteBaseLength());
		assertEquals(l.getLocalMaxBaseQual(), r.getRemoteMaxBaseQual());
		assertEquals(l.getLocalTotalBaseQual(), r.getRemoteTotalBaseQual());
		
		assertEquals(r.getLocalMapq(), l.getRemoteMapq());
		assertEquals(r.getLocalBaseLength(), l.getRemoteBaseLength());
		assertEquals(r.getLocalMaxBaseQual(), l.getRemoteMaxBaseQual());
		assertEquals(r.getLocalTotalBaseQual(), l.getRemoteTotalBaseQual());
	}
}
