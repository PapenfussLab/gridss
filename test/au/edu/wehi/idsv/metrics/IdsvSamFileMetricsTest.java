package au.edu.wehi.idsv.metrics;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import picard.analysis.InsertSizeMetrics;
import au.edu.wehi.idsv.TestHelper;

public class IdsvSamFileMetricsTest extends TestHelper {
	@Test
	public void shouldUseMADforStdDev() {
		IdsvMetrics im = new IdsvMetrics();
		InsertSizeMetrics ism = new InsertSizeMetrics();
		ism.MEDIAN_ABSOLUTE_DEVIATION = 10;
		RelevantMetrics metrics = new IdsvSamFileMetrics(ism, im);
		assertEquals(10 * 1.4826, metrics.getFragmentSizeStdDev(), 0.0001);
	}
	@Test
	public void shouldUseMaxProperPairSizeForMaxFragmentSize() {
		IdsvMetrics im = new IdsvMetrics();
		im.MAX_PROPER_PAIR_FRAGMENT_LENGTH = 10;
		InsertSizeMetrics ism = new InsertSizeMetrics();
		RelevantMetrics metrics = new IdsvSamFileMetrics(ism, im);
		assertEquals(10, metrics.getMaxFragmentSize());
	}
	@Test
	public void should_use_read_length_as_fragment_size_for_unpaired_reads() {
		IdsvMetrics im = new IdsvMetrics();
		im.MAX_PROPER_PAIR_FRAGMENT_LENGTH = 10;
		im.MAX_READ_LENGTH = 50;
		InsertSizeMetrics ism = new InsertSizeMetrics();
		IdsvSamFileMetrics metrics = new IdsvSamFileMetrics(ism, im);
		assertEquals(50, metrics.getMaxFragmentSize());
	}
}
