package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;

import java.util.List;

import htsjdk.samtools.SAMRecord;

import org.junit.Test;


public class SplitReadEvidenceTest extends TestHelper {
	@Test
	public void should_use_SA_tag_for_split_read() {
		SAMRecord r = Read(2, 1, "3S2M5S");
		//                    rrXmm-----
		r.setReadBases(    B("ACCTAGAGGG"));
		r.setBaseQualities(B("1234567890"));
		r.setMappingQuality(40);
		r.setAttribute("SA", "polyA,100,+,2M8S,10,0");
		SplitReadEvidence e = SplitReadEvidence.create(SES(), r).get(0);
		assertEquals(new BreakpointSummary(2, BWD, 1, 1, 0, FWD, 101, 101), e.getBreakendSummary());
		assertEquals("ACC", S(e.getBreakendSequence()));
		assertEquals("123", S(e.getBreakendQuality()));
		assertEquals("TA", S(e.getAnchorSequence()));
		assertEquals("45", S(e.getAnchorQuality()));
		assertEquals("C", e.getUntemplatedSequence());
		assertEquals(40, e.getLocalMapq());
		assertEquals(10, e.getRemoteMapq());
	}
	@Test
	public void should_return_adjacent_alignments() {
		SAMRecord r = Read(2, 1, "3S2M5S");
		//                    AB-mm-CDDD
		r.setReadBases(    B("ACCTAGAGGG"));
		r.setBaseQualities(B("1234567890"));
		r.setMappingQuality(40);
		r.setAttribute("SA", "polyA,300,+,6S1M3S,30,0;polyA,200,+,1S1M8S,20,0;polyA,400,+,7S3M,40,0;polyA,100,+,1M9S,10,0;");
		//                                C                  B                     D                  A
		List<SplitReadEvidence> list = SplitReadEvidence.create(SES(), r);
		assertEquals(2, list.size());
		assertEquals(new BreakpointSummary(2, BWD, 1, 1, 0, FWD, 200, 200), list.get(0).getBreakendSummary());
		assertEquals("TA", S(list.get(0).getAnchorSequence()));
		assertEquals("ACC", S(list.get(0).getBreakendSequence()));
		assertEquals("C", list.get(0).getUntemplatedSequence());
		assertEquals(new BreakpointSummary(2, FWD, 2, 2, 0, BWD, 300, 300), list.get(1).getBreakendSummary());
		assertEquals("TA", S(list.get(1).getAnchorSequence()));
		assertEquals("GAGGG", S(list.get(1).getBreakendSequence()));
		assertEquals("G", list.get(1).getUntemplatedSequence());
	}
	@Test
	public void should_handle_negative_strand_remote() {
		SAMRecord r = Read(2, 1, "3S2M5S");
		//                    -BBmm-CC--
		r.setReadBases(    B("ACCTAGAGGG"));
		r.setBaseQualities(B("1234567890"));
		r.setMappingQuality(40);
		r.setAttribute("SA", "polyA,100,-,7S2M1S,0,0;polyA,200,-,2S2M6S,0,0");
		List<SplitReadEvidence> list = SplitReadEvidence.create(SES(), r);
		assertEquals(new BreakpointSummary(2, BWD, 1, 1, 0, BWD, 100, 100), list.get(0).getBreakendSummary());
		assertEquals("TA", S(list.get(0).getAnchorSequence()));
		assertEquals("ACC", S(list.get(0).getBreakendSequence()));
		assertEquals("", list.get(0).getUntemplatedSequence());
		assertEquals(new BreakpointSummary(2, FWD, 2, 2, 0, FWD, 201, 201), list.get(1).getBreakendSummary());
		assertEquals("TA", S(list.get(1).getAnchorSequence()));
		assertEquals("GAGGG", S(list.get(1).getBreakendSequence()));
		assertEquals("G", list.get(1).getUntemplatedSequence());
	}
	@Test
	public void should_handle_negative_strand_local() {
		SAMRecord r = Read(2, 1, "5S2M3S");
		//                    --CC-mmBB-
		r.setReadBases(    B("ACCTAGAGGG"));
		r.setBaseQualities(B("1234567890"));
		r.setMappingQuality(40);
		r.setReadNegativeStrandFlag(true);
		r.setAttribute("SA", "polyA,100,+,1S2M7S,0,0;polyA,200,+,6S2M2S,0,0");
		List<SplitReadEvidence> list = SplitReadEvidence.create(SES(), r);
		assertEquals(new BreakpointSummary(2, FWD, 2, 2, 0, FWD, 101, 101), list.get(0).getBreakendSummary());
		assertEquals("GA", S(list.get(0).getAnchorSequence()));
		assertEquals("GGG", S(list.get(0).getBreakendSequence()));
		assertEquals("", list.get(0).getUntemplatedSequence());
		assertEquals(new BreakpointSummary(2, BWD, 1, 1, 0, BWD, 200, 200), list.get(1).getBreakendSummary());
		assertEquals("GA", S(list.get(1).getAnchorSequence()));
		assertEquals("ACCTA", S(list.get(1).getBreakendSequence()));
		assertEquals("A", list.get(1).getUntemplatedSequence());
	}
}
