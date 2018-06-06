package au.edu.wehi.idsv.alignment;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import org.junit.Test;
import org.junit.experimental.categories.Category;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.fastq.FastqRecord;

public class ExternalProcessStreamingAlignerTest {
	@Test
	@Category(ExternalAlignerTests.class)
	public void basic_pipes_test() throws IOException, InterruptedException {
		int COUNT = 1024;
		ExternalProcessStreamingAligner aligner = new ExternalProcessStreamingAligner(SamReaderFactory.makeDefault(), ExternalAlignerTests.COMMAND_LINE, ExternalAlignerTests.REFERENCE, 4);
		for (int i = 0; i < COUNT; i++) {
			aligner.asyncAlign(new FastqRecord(
					Integer.toString(i),
					"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
					Integer.toString(i),
					"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"));
		}
		aligner.close();
		for (int i = 0; i < COUNT; i++) {
			SAMRecord alignment = aligner.getAlignment();
			assertEquals(Integer.toString(i), alignment.getReadName());
		}
		aligner.close();
	}
	@Test
	@Category(ExternalAlignerTests.class)
	public void ping_pong_pipes_test() throws IOException, InterruptedException {
		int COUNT = 8;
		ExternalProcessStreamingAligner aligner = new ExternalProcessStreamingAligner(SamReaderFactory.makeDefault(), ExternalAlignerTests.COMMAND_LINE, ExternalAlignerTests.REFERENCE, 4);
		for (int i = 0; i < COUNT; i++) {
			aligner.asyncAlign(new FastqRecord(
					Integer.toString(i),
					"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
					Integer.toString(i),
					"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"));
			aligner.flush();
			SAMRecord alignment = aligner.getAlignment();
			assertEquals(Integer.toString(i), alignment.getReadName());
		}
		aligner.close();
	}
}
