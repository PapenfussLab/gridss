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
		assertEquals(new File("test.bam.gridss.working/test.bam.insert_size_metrics").getAbsolutePath(), C().getInsertSizeMetrics(new File("test.bam")).getAbsolutePath());
	}
	@Test
	public void getIdsvMetrics_should_use_empty_suffix() {
		assertEquals(new File("test.bam.gridss.working/test.bam.idsv_metrics").getAbsolutePath(), C().getIdsvMetrics(new File("test.bam")).getAbsolutePath());
	}
	@Test
	public void get_should_treat_working_files_as_belonging_to_their_parent() {
		assertEquals(new File("test.bam.gridss.working/test.bam.idsv_metrics").getAbsolutePath(), C().getIdsvMetrics(new File("test.bam.gridss.working/test.bamsv.bam")).getAbsolutePath());
	}
	private void testFileAssertMatch(String expected, File result) {
		assertEquals(new File("test.bam.gridss.working/" + expected).getAbsolutePath(), result.getAbsolutePath());
	}
	private static final File TEST_BAM = new File("test.bam");
	@Test
	public void should_match_constant() {
		testFileAssertMatch("test.bam.idsv_metrics", C().getIdsvMetrics(TEST_BAM));
		testFileAssertMatch("test.bam.insert_size_metrics", C().getInsertSizeMetrics(TEST_BAM));
		testFileAssertMatch("test.bam.realign.0.fq", C().getRealignmentFastq(TEST_BAM, 0));
		testFileAssertMatch("test.bam.realign.0.bam", C().getRealignmentBam(TEST_BAM, 0));
		testFileAssertMatch("test.bam.breakpoint.vcf", C().getBreakpointVcf(TEST_BAM));
	}
	@Test
	public void should_use_working_directory_if_set() throws IOException {
		File working = testFolder.newFolder("workingdir");
		File f = testFolder.newFile("test.bam");
		FileSystemContext fsc = new FileSystemContext(testFolder.getRoot(), working, 1);
		assertTrue(fsc.getIdsvMetrics(f).toString().contains("workingdir"));
	}
}
