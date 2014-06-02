package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMRecord;

import java.util.List;

import org.junit.Test;

public class ExtractEvidenceTest extends CommandLineTest {
	@Test
	public void should_generate_svbam_per_chr() {
		createInput();
		extractEvidence();
		workingFileShouldExist(".idsv.polyA.sv.bam");
		workingFileShouldExist(".idsv.polyACGT.sv.bam");
		workingFileShouldExist(".idsv.random.sv.bam");
		workingFileShouldExist(".idsv.Npower2.sv.bam");
		workingFileShouldExist(".idsv.REF.sv.bam");
		workingFileShouldExist(".idsv.REF2.sv.bam");
	}
	@Test
	public void sc_should_be_located_in_sv_bam_for_chr() {
		createInput(Read(1, 1, "50M50S"));
		extractEvidence();
		assertEquals(1, getRecords(".idsv.polyACGT.sv.bam").size());
	}
	@Test
	public void oea_should_be_located_in_sv_bam_for_chr() {
		createInput(OEA(1, 1, "100M", true));
		extractEvidence();
		assertEquals(1,  getRecords(".idsv.polyACGT.sv.bam").size());
	}
	@Test
	public void dp_should_be_located_in_sv_bam_for_chr() {
		createInput(DP(1, 1, "100M", true, 2, 2, "100M", true));
		extractEvidence();
		assertEquals(1,  getRecords(".idsv.polyACGT.sv.bam").size());
		assertEquals(1,  getRecords(".idsv.random.sv.bam").size());
	}
	@Test
	public void concordant_read_should_not_be_located_in_sv_bam() {
		createInput(Read(1, 1, "100M"));
		extractEvidence();
		assertEquals(0, getRecords(".idsv.polyACGT.sv.bam").size());
	}
	@Test
	public void oea_mate_should_be_located_in_sv_mate_bam_for_chr() {
		createInput(OEA(1, 1, "100M", true));
		extractEvidence();
		List<SAMRecord> rs = getRecords(".idsv.polyACGT.svmate.bam");
		assertEquals(1, rs.size());
		SAMRecord mate = rs.get(0);
		assertTrue(mate.getReadUnmappedFlag());
		assertEquals(1, (int)mate.getMateReferenceIndex());
		assertEquals(1, mate.getMateAlignmentStart());
	}
	@Test
	public void dp_mate_should_be_located_in_sv_mate_bam_for_chr() {
		createInput(DP(1, 1, "100M", true, 2, 2, "100M", true));
		extractEvidence();
		assertEquals(1,  getRecords(".idsv.polyACGT.svmate.bam").size());
		assertEquals(1,  getRecords(".idsv.random.svmate.bam").size());
		SAMRecord mate = getRecords(".idsv.polyACGT.sv.bam").get(0);
		assertFalse(mate.getReadUnmappedFlag());
		assertEquals(2, mate.getMateAlignmentStart());
		assertEquals(2, mate.getMateAlignmentStart());
	}
	@Test
	public void sv_mate_should_be_sorted_by_mate_coordinate() {
		createInput(DP(1, 1, "100M", true, 2, 5, "100M", true),
		   DP(1, 2, "100M", true, 2, 4, "100M", true),
		   DP(1, 3, "100M", true, 2, 6, "100M", true),
		   OEA(1, 4, "100M", false));
		extractEvidence();
		List<SAMRecord> rs = getRecords(".idsv.polyACGT.svmate.bam");
		assertEquals(4, rs.size());
		assertEquals(1, rs.get(0).getMateAlignmentStart());
		assertEquals(2, rs.get(1).getMateAlignmentStart());
		assertEquals(3, rs.get(2).getMateAlignmentStart());
		assertEquals(4, rs.get(3).getMateAlignmentStart());
		
		rs = getRecords(".idsv.random.svmate.bam"); 
		assertEquals(3, rs.size());
		assertEquals(4, rs.get(0).getMateAlignmentStart());
		assertEquals(5, rs.get(1).getMateAlignmentStart());
		assertEquals(6, rs.get(2).getMateAlignmentStart());
	}
	@Test
	public void should_create_metrics() {
		createInput();
		extractEvidence();
		workingFileShouldExist(".idsv.metrics.txt");
	}
	@Test
	public void should_set_NM_tag() {
		createInput(withSequence("TAAC", Read(0, 1, "1S3M")));
		extractEvidence();
		assertEquals(1, (int)getRecords(".idsv.polyA.sv.bam").get(0).getIntegerAttribute("NM"));
	}
}
