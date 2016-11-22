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
	public void getInsertSizeMetrics_should_use_empty_suffix() {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.metrics.insertsize").getAbsolutePath(), C().getInsertSizeMetrics(new File("test.bam")).getAbsolutePath());
	}
	@Test
	public void getIdsvMetrics_should_use_suffix() {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.metrics.idsv").getAbsolutePath(), C().getIdsvMetrics(new File("test.bam")).getAbsolutePath());
	}
	@Test
	public void get_should_treat_working_files_as_belonging_to_their_parent() {
		assertEquals(new File("test.bam.idsv.working/test.bam.idsv.metrics.idsv").getAbsolutePath(), C().getIdsvMetrics(new File("test.bam.idsv.working/test.bam.idsv.sv.bam")).getAbsolutePath());
	}
	private void testFileAssertMatch(String expected, File result) {
		assertEquals(new File("test.bam.idsv.working/" + expected).getAbsolutePath(), result.getAbsolutePath());
	}
	private static final File TEST_BAM = new File("test.bam");
	@Test
	public void should_match_constant() {
		testFileAssertMatch("test.bam.idsv.metrics.idsv", C().getIdsvMetrics(TEST_BAM));
		testFileAssertMatch("test.bam.idsv.metrics.insertsize", C().getInsertSizeMetrics(TEST_BAM));
		testFileAssertMatch("test.bam.idsv.realign.0.fq", C().getRealignmentFastq(TEST_BAM, 0));
		testFileAssertMatch("test.bam.idsv.realign.0.bam", C().getRealignmentBam(TEST_BAM, 0));
		testFileAssertMatch("test.bam.idsv.breakend.bam", C().getAssemblyRawBam(TEST_BAM));
		testFileAssertMatch("test.bam.idsv.breakpoint.vcf", C().getBreakpointVcf(TEST_BAM));
	}
	@Test
	public void test_remote() {
		testFileAssertMatch("test.bam.idsv.assembly.bam", C().getAssembly(TEST_BAM));
		testFileAssertMatch("test.bam.idsv.assemblymate.bam", C().getAssemblyMate(TEST_BAM));
		testFileAssertMatch("test.bam.idsv.scremote.bam", C().getSoftClipRemoteBam(TEST_BAM));
		testFileAssertMatch("test.bam.idsv.scremote.unsorted.bam", C().getSoftClipRemoteUnsortedBam(TEST_BAM));
		testFileAssertMatch("test.bam.idsv.realignremote.bam", C().getRealignmentRemoteBam(TEST_BAM));
		testFileAssertMatch("test.bam.idsv.realignremote.unsorted.bam", C().getRealignmentRemoteUnsortedBam(TEST_BAM));
	}
	@Test
	public void should_use_working_directory_if_set() throws IOException {
		File working = testFolder.newFolder("workingdir");
		File f = testFolder.newFile("test.bam");
		FileSystemContext fsc = new FileSystemContext(testFolder.getRoot(), working, 1);
		assertTrue(fsc.getAssembly(f).toString().contains("workingdir"));
	}
}
