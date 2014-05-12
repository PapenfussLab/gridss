package au.edu.wehi.socrates;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;

import org.junit.Test;


public class FileNamingConventionTest {
	@Test
	public void getMetrics_should_use_suffix() throws IOException {
		assertEquals(new File("test.bam.socrates.metrics.txt").getAbsolutePath(), FileNamingConvention.getMetrics(new File("test.bam")).getAbsolutePath());
	}
	@Test
	public void getMetrics_should_return_common_metrics_for_all_associated_files() throws IOException {
		assertEquals(new File("test.bam.socrates.metrics.txt").getAbsolutePath(), FileNamingConvention.getMetrics(new File("test.bam.socrates.sv.bam")).getAbsolutePath());
	}
	@Test
	public void getSVBam() throws IOException {
		assertEquals(new File("test.bam.socrates.sv.bam").getAbsolutePath(), FileNamingConvention.getSVBam(new File("test.bam")).getAbsolutePath());
	}
	@Test
	public void getSVBamForChr() throws IOException {
		assertEquals(new File("test.bam.socrates.chr1.sv.bam").getAbsolutePath(), FileNamingConvention.getSVBamForChr(new File("test.bam"), "chr1").getAbsolutePath());
	}
	@Test
	public void getMateBam() throws IOException {
		assertEquals(new File("test.bam.socrates.svmate.bam").getAbsolutePath(), FileNamingConvention.getMateBam(new File("test.bam")).getAbsolutePath());
	}
	@Test
	public void getMateBamForChr() throws IOException {
		assertEquals(new File("test.bam.socrates.chr1.svmate.bam").getAbsolutePath(), FileNamingConvention.getMateBamForChr(new File("test.bam"), "chr1").getAbsolutePath());
	}
	@Test
	public void getBreakendVcf() throws IOException {
		assertEquals(new File("test.bam.socrates.chr1.breakend.vcf").getAbsolutePath(), FileNamingConvention.getBreakendVcfForChr(new File("test.bam"), "chr1").getAbsolutePath());
	}
	@Test
	public void getBreakpointVcf() throws IOException {
		assertEquals(new File("test.bam.socrates.chr1-chr2.breakpoint.vcf").getAbsolutePath(), FileNamingConvention.getBreakpointVcf(new File("test.bam"), "chr1", "chr2").getAbsolutePath());
	}
	@Test
	public void getOutputVcf() throws IOException {
		assertEquals(new File("test.bam.socrates.vcf").getAbsolutePath(), FileNamingConvention.getOutputVcf(new File("test.bam")).getAbsolutePath());
	}
}
