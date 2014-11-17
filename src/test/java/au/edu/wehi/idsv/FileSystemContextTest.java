package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

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
	public void getInsertSizeMetrics_should_use_suffix() {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.metrics.insertsize.txt").getAbsolutePath(), C().getInsertSizeMetrics(new File("test.bam")).getAbsolutePath());
	}
	@Test
	public void getIdsvMetrics_should_use_suffix() {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.metrics.idsv.txt").getAbsolutePath(), C().getIdsvMetrics(new File("test.bam")).getAbsolutePath());
	}
	@Test
	public void get_should_treat_working_files_as_belonging_to_their_parent() {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.metrics.idsv.txt").getAbsolutePath(), C().getIdsvMetrics(new File("test.bam.idsv.working/test.bam.idsv.sv.bam")).getAbsolutePath());
	}
	@Test
	public void getReadPairBam() {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.rp.bam").getAbsolutePath(), C().getReadPairBam(new File("test.bam")).getAbsolutePath());
	}
	@Test
	public void getReadPairBamForChr() {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.chr1.rp.bam").getAbsolutePath(), C().getReadPairBamForChr(new File("test.bam"), "chr1").getAbsolutePath());
	}
	@Test
	public void getSoftClipBam() {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.sc.bam").getAbsolutePath(), C().getSoftClipBam(new File("test.bam")).getAbsolutePath());
	}
	@Test
	public void getSoftClipBamForChr() {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.chr1.sc.bam").getAbsolutePath(), C().getSoftClipBamForChr(new File("test.bam"), "chr1").getAbsolutePath());
	}
	@Test
	public void getMateBam() {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.rpmate.bam").getAbsolutePath(), C().getMateBam(new File("test.bam")).getAbsolutePath());
	}
	@Test
	public void getMateBamForChr() {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.chr1.rpmate.bam").getAbsolutePath(), C().getMateBamForChr(new File("test.bam"), "chr1").getAbsolutePath());
	}
	@Test
	public void getBreakendVcf() {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.chr1.breakend.bam").getAbsolutePath(), C().getAssemblyRawBamForChr(new File("test.bam"), "chr1").getAbsolutePath());
	}
	@Test
	public void getBreakpointVcf() {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.chr1-chr2.breakpoint.vcf").getAbsolutePath(), C().getBreakpointVcf(new File("test.bam"), "chr1", "chr2").getAbsolutePath());
	}
	private void testFileAssertMatch(String expected, File result) {
		assertEquals(new File("test.bam.idsv.working/" + expected).getAbsolutePath(), result.getAbsolutePath());
	}
	private static final File TEST_BAM = new File("test.bam");
	@Test
	public void should_match_constant() {
		testFileAssertMatch("test.bam.idsv.rp.bam", C().getReadPairBam(TEST_BAM));
		testFileAssertMatch("test.bam.idsv.sc.bam", C().getSoftClipBam(TEST_BAM));
		testFileAssertMatch("test.bam.idsv.rpmate.bam", C().getMateBam(TEST_BAM));
		testFileAssertMatch("test.bam.idsv.metrics.idsv.txt", C().getIdsvMetrics(TEST_BAM));
		testFileAssertMatch("test.bam.idsv.metrics.insertsize.txt", C().getInsertSizeMetrics(TEST_BAM));
		testFileAssertMatch("test.bam.idsv.realign.fq", C().getRealignmentFastq(TEST_BAM));
		testFileAssertMatch("test.bam.idsv.realign.bam", C().getRealignmentBam(TEST_BAM));
		testFileAssertMatch("test.bam.idsv.breakend.bam", C().getAssemblyRawBam(TEST_BAM));
		testFileAssertMatch("test.bam.idsv.breakpoint.vcf", C().getBreakpointVcf(TEST_BAM));
	}
	@Test
	public void should_match_constant_per_chr() {
		testFileAssertMatch("test.bam.idsv.chr.rp.bam", C().getReadPairBamForChr(TEST_BAM, "chr"));
		testFileAssertMatch("test.bam.idsv.chr.sc.bam", C().getSoftClipBamForChr(TEST_BAM, "chr"));
		testFileAssertMatch("test.bam.idsv.chr.rpmate.bam", C().getMateBamForChr(TEST_BAM, "chr"));
		testFileAssertMatch("test.bam.idsv.chr.realign.fq", C().getRealignmentFastqForChr(TEST_BAM, "chr"));
		testFileAssertMatch("test.bam.idsv.chr.realign.bam", C().getRealignmentBamForChr(TEST_BAM, "chr"));
		testFileAssertMatch("test.bam.idsv.chr.breakend.bam", C().getAssemblyRawBamForChr(TEST_BAM, "chr"));
		testFileAssertMatch("test.bam.idsv.chr1-chr2.breakpoint.vcf", C().getBreakpointVcf(TEST_BAM, "chr1", "chr2"));
	}
	@Test
	public void test_remote() {
		testFileAssertMatch("test.bam.idsv.assemblyremote.bam", C().getAssembly(TEST_BAM));
		testFileAssertMatch("test.bam.idsv.assemblyremote.unsorted.bam", C().getAssemblyMate(TEST_BAM));
		testFileAssertMatch("test.bam.idsv.scremote.bam", C().getSoftClipRemoteBam(TEST_BAM));
		testFileAssertMatch("test.bam.idsv.scremote.unsorted.bam", C().getSoftClipRemoteUnsortedBam(TEST_BAM));
		testFileAssertMatch("test.bam.idsv.realignremote.bam", C().getRealignmentRemoteBam(TEST_BAM));
		testFileAssertMatch("test.bam.idsv.realignremote.unsorted.bam", C().getRealignmentRemoteUnsortedBam(TEST_BAM));
	}
	@Test
	public void test_remote_per_chr() {
		testFileAssertMatch("test.bam.idsv.chr.assemblyremote.bam", C().getAssemblyForChr(TEST_BAM, "chr"));
		testFileAssertMatch("test.bam.idsv.chr.assemblyremote.unsorted.bam", C().getAssemblyMateForChr(TEST_BAM, "chr"));
		testFileAssertMatch("test.bam.idsv.chr.scremote.bam", C().getSoftClipRemoteBamForChr(TEST_BAM, "chr"));
		testFileAssertMatch("test.bam.idsv.chr.scremote.unsorted.bam", C().getSoftClipRemoteUnsortedBamForChr(TEST_BAM, "chr"));
		testFileAssertMatch("test.bam.idsv.chr.realignremote.bam", C().getRealignmentRemoteBamForChr(TEST_BAM, "chr"));
		testFileAssertMatch("test.bam.idsv.chr.realignremote.unsorted.bam", C().getRealignmentRemoteUnsortedBamForChr(TEST_BAM, "chr"));
	}
	@Test
	public void should_use_working_directory_if_set() throws IOException {
		File working = testFolder.newFolder("working");
		File f = testFolder.newFile("test.bam");
		FileSystemContext fsc = new FileSystemContext(testFolder.getRoot(), working, 1);
		assertTrue(fsc.getReadPairBam(f).toString().contains("working"));
	}
}
