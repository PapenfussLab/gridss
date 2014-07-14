package au.edu.wehi.idsv;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;


public class FileSystemContextTest {
	@Rule
    public TemporaryFolder testFolder = new TemporaryFolder();
	public FileSystemContext C() {
		return new FileSystemContext(testFolder.getRoot(), 1);
	}
	@Test
	public void getMetrics_should_use_suffix() {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.metrics.txt").getAbsolutePath(), C().getMetrics(new File("test.bam")).getAbsolutePath());
	}
	@Test
	public void get_should_treat_working_files_as_belonging_to_their_parent() {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.metrics.txt").getAbsolutePath(), C().getMetrics(new File("test.bam.idsv.working/test.bam.idsv.sv.bam")).getAbsolutePath());
	}
	@Test
	public void getSVBam() {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.sv.bam").getAbsolutePath(), C().getSVBam(new File("test.bam")).getAbsolutePath());
	}
	@Test
	public void getSVBamForChr() {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.chr1.sv.bam").getAbsolutePath(), C().getSVBamForChr(new File("test.bam"), "chr1").getAbsolutePath());
	}
	@Test
	public void getMateBam() {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.svmate.bam").getAbsolutePath(), C().getMateBam(new File("test.bam")).getAbsolutePath());
	}
	@Test
	public void getMateBamForChr() {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.chr1.svmate.bam").getAbsolutePath(), C().getMateBamForChr(new File("test.bam"), "chr1").getAbsolutePath());
	}
	@Test
	public void getBreakendVcf() {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.chr1.breakend.vcf").getAbsolutePath(), C().getBreakendVcfForChr(new File("test.bam"), "chr1").getAbsolutePath());
	}
	@Test
	public void getBreakpointVcf() {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.chr1-chr2.breakpoint.vcf").getAbsolutePath(), C().getBreakpointlVcf(new File("test.bam"), "chr1", "chr2").getAbsolutePath());
	}
	private void testbamAssertMatch(String expected, File result) {
		assertEquals(new File("test.bam.idsv.working/" + expected).getAbsolutePath(), result.getAbsolutePath());
	}
	private static final File TEST_BAM = new File("test.bam");
	@Test
	public void should_match_constant() {
		testbamAssertMatch("test.bam.idsv.sv.bam", C().getSVBam(TEST_BAM));
		testbamAssertMatch("test.bam.idsv.svmate.bam", C().getMateBam(TEST_BAM));
		testbamAssertMatch("test.bam.idsv.metrics.txt", C().getMetrics(TEST_BAM));
		testbamAssertMatch("test.bam.idsv.realign.fq", C().getRealignmentFastq(TEST_BAM));
		testbamAssertMatch("test.bam.idsv.realign.bam", C().getRealignmentBam(TEST_BAM));
		testbamAssertMatch("test.bam.idsv.breakend.vcf", C().getBreakendVcf(TEST_BAM));
		testbamAssertMatch("test.bam.idsv.breakpoint.vcf", C().getBreakpointVcf(TEST_BAM));
	}
	@Test
	public void should_match_constant_per_chr() {
		testbamAssertMatch("test.bam.idsv.chr.sv.bam", C().getSVBamForChr(TEST_BAM, "chr"));
		testbamAssertMatch("test.bam.idsv.chr.svmate.bam", C().getMateBamForChr(TEST_BAM, "chr"));
		testbamAssertMatch("test.bam.idsv.chr.realign.fq", C().getRealignmentFastqForChr(TEST_BAM, "chr"));
		testbamAssertMatch("test.bam.idsv.chr.realign.bam", C().getRealignmentBamForChr(TEST_BAM, "chr"));
		testbamAssertMatch("test.bam.idsv.chr.breakend.vcf", C().getBreakendVcfForChr(TEST_BAM, "chr"));
		testbamAssertMatch("test.bam.idsv.chr1-chr2.breakpoint.vcf", C().getBreakpointlVcf(TEST_BAM, "chr1", "chr2"));
	}
	@Test
	public void should_use_working_directory_if_set() throws IOException {
		File working = testFolder.newFolder("working");
		File f = testFolder.newFile("test.bam");
		FileSystemContext fsc = new FileSystemContext(testFolder.getRoot(), working, 1);
		assertTrue(fsc.getSVBam(f).toString().contains("working"));
	}
}
