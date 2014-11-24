package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;

import org.junit.Test;

import com.google.common.io.Files;


public class EvidenceSourceTest extends IntermediateFilesTest {
	public class TestEvidenceSource extends EvidenceSource {
		public TestEvidenceSource(boolean perChr, File input) {
			super(getCommandlineContext(perChr), input);
		}
		protected boolean processed = false;
		@Override
		public int getMaxConcordantFragmentSize() {
			return 300;
		}
	}
	@Test
	public void getFileIntermediateDirectoryBasedOn_should_be_input_file() {
		EvidenceSource es = new TestEvidenceSource(false, input);
		assertEquals(input, es.getFileIntermediateDirectoryBasedOn());
	}
	@Test
	public void getRealignmentScript_should_not_write_realignment_for_missing_fastq() {
		for (boolean perChr : new boolean[] {true, false} ) {
			EvidenceSource es = new TestEvidenceSource(perChr, input);
			assertEquals("", es.getRealignmentScript(1));
		}
	}
	@Test
	public void getRealignmentScript_should_not_write_realignment_for_existing_bam() throws IOException {
		FileSystemContext fsc = getCommandlineContext(false).getFileSystemContext();
		EvidenceSource es = new TestEvidenceSource(false, input);
		
		Files.write(new byte[] { '>' }, fsc.getRealignmentFastq(input));
		Files.write(new byte[] { '#' }, fsc.getRealignmentBam(input));
		assertEquals("", es.getRealignmentScript(1));
	}
	@Test
	public void getRealignmentScript_should_not_write_realignment_empty_fastq() throws IOException {
		FileSystemContext fsc = getCommandlineContext(false).getFileSystemContext();
		EvidenceSource es = new TestEvidenceSource(false, input);		
		Files.write(new byte[0], fsc.getRealignmentFastq(input));
		assertEquals("", es.getRealignmentScript(1));
	}
	@Test
	public void getRealignmentScript_should_write_realignment_for_fastq_records() throws IOException {
		FileSystemContext fsc = getCommandlineContext(false).getFileSystemContext();
		EvidenceSource es = new TestEvidenceSource(false, input);		
		Files.write(new byte[] { '>' }, fsc.getRealignmentFastq(input));
		assertNotEquals("", es.getRealignmentScript(1));
	}
	@Test
	public void getRealignmentScript_should_write_realignment_for_out_of_date_bam() throws IOException, InterruptedException {
		FileSystemContext fsc = getCommandlineContext(false).getFileSystemContext();
		EvidenceSource es = new TestEvidenceSource(false, input);
		Files.write(new byte[] { '#' }, fsc.getRealignmentBam(input));
		Thread.sleep(1); // wait 1ms so timestamps are different
		Files.write(new byte[] { '>' }, fsc.getRealignmentFastq(input));
		assertNotEquals("", es.getRealignmentScript(1));
	}
	@Test
	public void getRealignmentScript_should_default_to_bowtie2() throws IOException {
		FileSystemContext fsc = getCommandlineContext(false).getFileSystemContext();
		EvidenceSource es = new TestEvidenceSource(false, input);
		Files.write(new byte[] { '>' }, fsc.getRealignmentFastq(input));
		String script = es.getRealignmentScript(1);
		assertTrue(script.contains("bowtie2"));
	}
}
