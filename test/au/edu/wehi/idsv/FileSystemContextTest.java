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
	public FileSystemContext C() throws IOException {
		return new FileSystemContext(testFolder.getRoot(), 1);
	}
	@Test
	public void getMetrics_should_use_suffix() throws IOException {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.metrics.txt").getAbsolutePath(), C().getMetrics(new File("test.bam")).getAbsolutePath());
	}
	@Test
	public void get_should_treat_working_files_as_belonging_to_their_parent() throws IOException {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.metrics.txt").getAbsolutePath(), C().getMetrics(new File("test.bam.idsv.working/test.bam.idsv.sv.bam")).getAbsolutePath());
	}
	@Test
	public void getSVBam() throws IOException {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.sv.bam").getAbsolutePath(), C().getSVBam(new File("test.bam")).getAbsolutePath());
	}
	@Test
	public void getSVBamForChr() throws IOException {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.chr1.sv.bam").getAbsolutePath(), C().getSVBamForChr(new File("test.bam"), "chr1").getAbsolutePath());
	}
	@Test
	public void getMateBam() throws IOException {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.svmate.bam").getAbsolutePath(), C().getMateBam(new File("test.bam")).getAbsolutePath());
	}
	@Test
	public void getMateBamForChr() throws IOException {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.chr1.svmate.bam").getAbsolutePath(), C().getMateBamForChr(new File("test.bam"), "chr1").getAbsolutePath());
	}
	@Test
	public void getBreakendVcf() throws IOException {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.chr1.breakend.vcf").getAbsolutePath(), C().getBreakendVcfForChr(new File("test.bam"), "chr1").getAbsolutePath());
	}
	@Test
	public void getBreakpointVcf() throws IOException {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.chr1-chr2.breakpoint.vcf").getAbsolutePath(), C().getRawCallVcf(new File("test.bam"), "chr1", "chr2").getAbsolutePath());
	}
	@Test
	public void should_use_working_directory_if_set() throws IOException {
		File working = testFolder.newFolder("working");
		File f = testFolder.newFile("test.bam");
		FileSystemContext fsc = new FileSystemContext(testFolder.getRoot(), working, 1);
		assertTrue(fsc.getSVBam(f).toString().contains("working"));
	}
}
