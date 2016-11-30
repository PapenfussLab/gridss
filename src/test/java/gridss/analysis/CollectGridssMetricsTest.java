package gridss.analysis;

import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;

import org.junit.Ignore;
import org.junit.Test;

import au.edu.wehi.idsv.IntermediateFilesTest;

public class CollectGridssMetricsTest extends IntermediateFilesTest {
	@Test
	public void should_generate_all_metrics() throws IOException {
		String prefix = new File(testFolder.getRoot(), "output").getAbsolutePath();
		CollectGridssMetrics collect = new CollectGridssMetrics();
		collect.instanceMain(new String[] {
			"INPUT=" + new File("src/test/resources/203541.bam").getAbsolutePath(),
			"OUTPUT=" + prefix,
		});
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
		CollectGridssMetrics collect = new CollectGridssMetrics();
		collect.instanceMain(new String[] {
			"INPUT=" + new File("src/test/resources/203541.bam").getAbsolutePath(),
			"OUTPUT=" + prefix,
		});
		assertTrue(new File(prefix + ".mapq_histogram.pdf").exists());
	}
}

