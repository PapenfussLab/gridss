package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMRecord;

import java.util.List;

import org.junit.Test;

public class IndelEvidenceTest extends TestHelper {
	@Test
	public void create_should_generate_both_sides() {
		SAMRecord r = Read(2, 1, "5M5D5M");
		r.setMappingQuality(40);
		List<IndelEvidence> e = IndelEvidence.create(SES(), r);
		assertEquals(2, e.size());
		assertEquals(e.get(0), e.get(1).asRemote());
		assertEquals(e.get(0).asRemote(), e.get(1));
	}
	@Test
	public void create_should_generate_all_indels() {
		SAMRecord r = Read(2, 1, "5M1D5M2I2D5M3I5M");
		r.setMappingQuality(40);
		List<IndelEvidence> e = IndelEvidence.create(SES(), r);
		assertEquals(6, e.size());
	}
	@Test
	public void qual_should_use_library_indel_distribution_and_split_score_to_both_sides() {
		SAMRecord r = Read(2, 1, "5M5D5M");
		r.setMappingQuality(40);
		IndelEvidence e = IndelEvidence.create(SES(), r).get(0);
		assertEquals(5, e.getBreakendQual() + e.asRemote().getBreakendQual(), 0.01);
		assertEquals(5, e.getBreakpointQual() + e.asRemote().getBreakpointQual(), 0.01);
	}
	@Test
	public void should_allow_source_read_alignment_to_either_strand() {
		SAMRecord r = Read(2, 1, "2M2D2M");
		r.setMappingQuality(40);
		r.setReadNegativeStrandFlag(true);
		IndelEvidence e = IndelEvidence.create(SES(), r).get(0);
		assertEquals(new BreakpointSummary(2, FWD, 2, 2, 2, BWD, 5, 5), e.getBreakendSummary());
		assertEquals(new BreakpointSummary(2, BWD, 5, 5, 2, FWD, 2, 2), e.asRemote().getBreakendSummary());
	}
	@Test
	public void should_span_deletion_event() {
		SAMRecord r = Read(2, 1, "2M2D2M");
		r.setReadBases(B("NNNN"));
		r.setMappingQuality(40);
		IndelEvidence e = IndelEvidence.create(SES(), r).get(0);
		assertEquals(new BreakpointSummary(2, FWD, 2, 2, 2, BWD, 5, 5), e.getBreakendSummary());
		assertEquals(new BreakpointSummary(2, BWD, 5, 5, 2, FWD, 2, 2), e.asRemote().getBreakendSummary());
	}
	@Test
	public void should_span_insertion_event() {
		SAMRecord r = Read(2, 1, "2M2I2M");
		r.setReadBases(B("NNTTNN"));
		r.setMappingQuality(40);
		IndelEvidence e = IndelEvidence.create(SES(), r).get(0);
		assertEquals(new BreakpointSummary(2, FWD, 2, 2, 2, BWD, 3, 3), e.getBreakendSummary());
		assertEquals(new BreakpointSummary(2, BWD, 3, 3, 2, FWD, 2, 2), e.asRemote().getBreakendSummary());
		assertEquals("TT", e.getUntemplatedSequence());
		assertEquals("TT", e.asRemote().getUntemplatedSequence());
	}
	@Test
	public void indels_should_have_unique_evidenceID() {
		SAMRecord r = Read(2, 1, "2M2I1M2D2M");
		r.setReadBases(B("NNNNNNN"));
		r.setMappingQuality(40);
		List<IndelEvidence> e = IndelEvidence.create(SES(), r);
		assertEquals(4, e.size());
		assertEquals(4, e.stream().map(ie -> ie.getEvidenceID()).distinct().count());
	}
	@Test
	public void getBreakendSummary_should_return_location() {
		// 1234567890
		// MMDDDMMMM
		SAMRecord r = Read(2, 1, "2M3D4M");
		List<IndelEvidence> e = IndelEvidence.create(SES(), r);
		assertEquals(new BreakpointSummary(2, FWD, 2, 2, 2, BWD, 6, 6), e.get(0).getBreakendSummary());
		assertEquals(new BreakpointSummary(2, BWD, 6, 6, 2, FWD, 2, 2), e.get(1).getBreakendSummary());
	}
	
	@Test
	public void getBreakendSequence_should_return_inserted_remote() {
		// 1234567890
		// MMDDDMMMM
		SAMRecord r = Read(2, 1, "2M3I1D1M");
		r.setReadBases(B("ACGTTA"));
		List<IndelEvidence> e = IndelEvidence.create(SES(), r);
		assertEquals("GTTA", S(e.get(0).getBreakendSequence()));
		assertEquals("ACGTT", S(e.get(1).getBreakendSequence()));
	}
	
	@Test
	public void getBreakendQuality_should_return_inserted_remote() {
		// 1234567890
		// MMDDDMMMM
		SAMRecord r = Read(2, 1, "2M3I1D1M");
		r.setReadBases(B("ACGTTA"));
		r.setBaseQualities(B("123456"));
		List<IndelEvidence> e = IndelEvidence.create(SES(), r);
		assertEquals("3456", S(e.get(0).getBreakendQuality()));
		assertEquals("12345", S(e.get(1).getBreakendQuality()));
	}
	
	@Test
	public void getAnchorSequence_should_return_local() {
		// 1234567890
		// MMDDDMMMM
		SAMRecord r = Read(2, 1, "2M3I1D1M");
		r.setReadBases(B("ACGTTA"));
		r.setBaseQualities(B("123456"));
		List<IndelEvidence> e = IndelEvidence.create(SES(), r);
		assertEquals("AC", S(e.get(0).getAnchorSequence()));
		assertEquals("A", S(e.get(1).getAnchorSequence()));
	}
	
	@Test
	public void getAnchorQuality_should_return_local() {
		// 1234567890
		// MMDDDMMMM
		SAMRecord r = Read(2, 1, "2M3I1D1M");
		r.setReadBases(B("ACGTTA"));
		r.setBaseQualities(B("123456"));
		List<IndelEvidence> e = IndelEvidence.create(SES(), r);
		assertEquals("12", S(e.get(0).getAnchorQuality()));
		assertEquals("6", S(e.get(1).getAnchorQuality()));
	}
	
	@Test
	public void getMapq_should_match_read_mapq() {
		SAMRecord r = Read(2, 1, "2M2I1M2D2M");
		r.setMappingQuality(40);
		List<IndelEvidence> e = IndelEvidence.create(SES(), r);
		assertEquals(40, e.get(0).getLocalMapq());
		assertEquals(40, e.get(0).getRemoteMapq());
	}
	
	@Test
	public void isBreakendExact() {
		SAMRecord r = Read(2, 1, "2M2I1M2D2M");
		r.setMappingQuality(40);
		List<IndelEvidence> e = IndelEvidence.create(SES(), r);
		assertTrue(e.get(0).isBreakendExact());
	}
	@Test
	public void getUntemplatedSequence_should_return_insert() {
		SAMRecord r = Read(2, 1, "2M3I1D1M");
		r.setReadBases(B("ACGTTA"));
		r.setBaseQualities(B("123456"));
		List<IndelEvidence> e = IndelEvidence.create(SES(), r);
		assertEquals("GTT", e.get(0).getUntemplatedSequence());
		assertEquals("GTT", e.get(1).getUntemplatedSequence());
	}
}