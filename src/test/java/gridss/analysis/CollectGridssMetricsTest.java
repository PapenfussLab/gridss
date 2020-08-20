package gridss.analysis;

import au.edu.wehi.idsv.IntermediateFilesTest;
import gridss.cmdline.CommandLineProgramHelper;
import org.junit.Ignore;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

import static org.junit.Assert.assertTrue;

public class CollectGridssMetricsTest extends IntermediateFilesTest {
	@Test
	public void should_generate_all_metrics() throws IOException {
		String prefix = new File(testFolder.getRoot(), "output").getAbsolutePath();
		CommandLineProgramHelper cmd = new CommandLineProgramHelper(new CollectGridssMetrics());
		cmd.addArg("I", new File("src/test/resources/203541.bam").getAbsolutePath());
		cmd.addArg("O", prefix);
		cmd.addArg("THRESHOLD_COVERAGE", 1000);
		cmd.run();
		assertTrue(new File(prefix + ".cigar_metrics").exists());
		assertTrue(new File(prefix + ".insert_size_metrics").exists());
		assertTrue(new File(prefix + ".mapq_metrics").exists());
		assertTrue(new File(prefix + ".idsv_metrics").exists());
		assertTrue(new File(prefix + ".tag_metrics").exists());
	}
	@Test
	@Ignore("Replaced Rscript with placeholder noop executable to reduce unit test runtime")
	public void should_generate_histogram() throws IOException {
		String prefix = new File(testFolder.getRoot(), "output").getAbsolutePath();
		CommandLineProgramHelper cmd = new CommandLineProgramHelper(new CollectGridssMetrics());
		cmd.addArg("I", new File("src/test/resources/203541.bam").getAbsolutePath());
		cmd.addArg("O", prefix);
		cmd.addArg("THRESHOLD_COVERAGE", 1000);
		cmd.run();
		assertTrue(new File(prefix + ".mapq_histogram.pdf").exists());
	}
}

