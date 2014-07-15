package au.edu.wehi.idsv.metrics;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;


public class IdsvMetricsCollectorTest extends TestHelper {
	@Test
	public void should_calc_max_read_length() {
		IdsvMetricsCollector c = new IdsvMetricsCollector();
		c.acceptRecord(RP(0, 1, 2, 1)[0], null);
		c.acceptRecord(RP(0, 1, 2, 1)[1], null);
		c.acceptRecord(RP(0, 1, 7, 5)[0], null);
		c.acceptRecord(RP(0, 1, 7, 5)[1], null);
		c.acceptRecord(Read(0, 1, "100M"), null);
		MetricsFile<IdsvMetrics, Integer> imcf = new MetricsFile<IdsvMetrics, Integer>();
		c.addAllLevelsToFile(imcf);
		assertEquals(100, (int)((IdsvMetrics)imcf.getMetrics().get(0)).MAX_READ_LENGTH);
	}
	@Test
	public void should_calc_max_fragment_size() {
		IdsvMetricsCollector c = new IdsvMetricsCollector();
		c.acceptRecord(Read(0, 1, "100M"), null);
		c.acceptRecord(RP(0, 1, 2, 1)[0], null);
		c.acceptRecord(RP(0, 1, 7, 5)[0], null);
		SAMRecord r = RP(0, 1, 100, 5)[0]; // should ignore this as it's not a proper pair
		r.setProperPairFlag(false);
		c.acceptRecord(r, null);
		
		// 12345678901234567890
		// ----> <----
		c.acceptRecord(Read(0, 1, "100M"), null);
		MetricsFile<IdsvMetrics, Integer> imcf = new MetricsFile<IdsvMetrics, Integer>();
		c.addAllLevelsToFile(imcf);
		assertEquals(11, (int)((IdsvMetrics)imcf.getMetrics().get(0)).MAX_PROPER_PAIR_FRAGMENT_LENGTH);
		
		c = new IdsvMetricsCollector();
		r = RP(0, 1, 100, 5)[1]; // mate before
		c.acceptRecord(r, null);
		c.addAllLevelsToFile(imcf);
		assertEquals(11, (int)((IdsvMetrics)imcf.getMetrics().get(0)).MAX_PROPER_PAIR_FRAGMENT_LENGTH);
	}
	@Test
	public void should_add_metrics_to_file() {
		IdsvMetricsCollector c = new IdsvMetricsCollector();
		MetricsFile<IdsvMetrics, Integer> imcf = new MetricsFile<IdsvMetrics, Integer>();
		c.addAllLevelsToFile(imcf);
		assertEquals(1, imcf.getMetrics().size());
		assertTrue(imcf.getMetrics().get(0) instanceof IdsvMetrics);
	}
}
