package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

import org.junit.Before;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;


public abstract class RemoteEvidenceTest extends TestHelper  {
	protected int anchorLength(final String allBases, final String realignCigar) {
		return allBases.length() - (new SAMRecord(null) {{ setCigarString(realignCigar); }}).getCigar().getReadLength();
	}
	protected int realignReferenceIndex;
	protected int realignAlignmentStart;
	@Before
	public void setup() {
		realignReferenceIndex = 1;
		realignAlignmentStart = 100;
	}
	protected SAMRecord realignSAM(final BreakendSummary bs, final String allBases, final String realignCigar, final boolean realignNegativeStrand) {
		return new SAMRecord(getContext().getBasicSamHeader()) {{
			setReferenceIndex(realignReferenceIndex);
			setAlignmentStart(realignAlignmentStart);
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
	@Test
	public void should_have_different_evidenceID() {
		assertNotEquals(L().getEvidenceID(), R().getEvidenceID());
	}
	@Test
	public void should_prefix_evidenceID_with_R() {
		assertEquals("R" + L().getEvidenceID(), R().getEvidenceID());
	}
	@Test
	public void should_prefix_toString_with_R() {
		assertEquals("R", R().toString().substring(0, 1));
	}
	@Test
	public void getBreakpointQual_should_be_unchanged() {
		assertEquals(L().getBreakpointQual(), R().getBreakpointQual(), 0);
	}
	@Test
	public void getBreakendQual_should_be_unchanged() {
		assertEquals(L().getBreakendQual(), R().getBreakendQual(), 0);
	}
	@Test
	public void local_remote_should_be_same_breakpoint() {
		assertEquals(L().getBreakendSummary(), R().getBreakendSummary().remoteBreakpoint());
	}
	@Test
	public void homology_sequence_should_be_as_if_local() {
		realignReferenceIndex = 0;
		DirectedBreakpoint bp = makeRemote(new BreakendSummary(0, FWD, 100, 100), "AAAAAAAA", "1M", false);
		assertEquals("AAAAAAAA", bp.getHomologySequence());
		assertEquals(1, bp.getHomologyAnchoredBaseCount());
		
		realignReferenceIndex = 0;
		bp = makeRemote(new BreakendSummary(0, BWD, 100, 100), "AAAAAAAA", "1M", false);
		assertEquals("AAAAAAAA", bp.getHomologySequence());
		assertEquals(1, bp.getHomologyAnchoredBaseCount());
		
		 // 3519-3523 = TTTTT
		realignReferenceIndex = 2;
		realignAlignmentStart = 3519;
		bp = makeRemote(new BreakendSummary(0, FWD, 100, 100), "AAAAA", "1M", true);
		assertEquals("TTTTT", bp.getHomologySequence());
		assertEquals(1, bp.getHomologyAnchoredBaseCount());
		
		realignReferenceIndex = 2;
		realignAlignmentStart = 3523;
		bp = makeRemote(new BreakendSummary(0, BWD, 100, 100), "AAAAA", "1M", true);
		assertEquals("TTTTT", bp.getHomologySequence());
		assertEquals(1, bp.getHomologyAnchoredBaseCount());
		
		realignReferenceIndex = 1;
		realignAlignmentStart = 3;
		bp = makeRemote(new BreakendSummary(1, FWD, 1, 1), "AC", "1M", true);
		assertEquals("GT", bp.getHomologySequence());
		assertEquals(1, bp.getHomologyAnchoredBaseCount());
	}
}
