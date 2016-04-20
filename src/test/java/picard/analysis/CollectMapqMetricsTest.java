package picard.analysis;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;

import org.junit.Test;

import au.edu.wehi.idsv.IntermediateFilesTest;


public class CollectMapqMetricsTest extends IntermediateFilesTest {
	@Test
	public void should_generate_metrics() throws IOException {
		CollectMapqMetrics cmm = new CollectMapqMetrics();
		cmm.INPUT = new File("src/test/resources/203541.bam");
		cmm.OUTPUT = new File(testFolder.getRoot(), "mapqmetrics.txt");
		cmm.Histogram_FILE = new File(testFolder.getRoot(), "mapqhistogram.pdf");
		cmm.doWork();
		assertTrue(cmm.OUTPUT.isFile());
	}
	@Test
	public void should_generate_histogram() throws IOException {
		CollectMapqMetrics cmm = new CollectMapqMetrics();
		File histogram = new File(testFolder.getRoot(), "mapqhistogram.pdf");
		cmm.INPUT = new File("src/test/resources/203541.bam");
		cmm.OUTPUT = new File(testFolder.getRoot(), "mapqmetrics.txt");
		cmm.Histogram_FILE = histogram;
		cmm.doWork();
		assertTrue(histogram.isFile());
	}
}
