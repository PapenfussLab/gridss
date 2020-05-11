package gridss.analysis;

import au.edu.wehi.idsv.IntermediateFilesTest;
import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.Histogram.Bin;
import org.junit.Ignore;
import org.junit.Test;
import picard.analysis.SinglePassSamProgram;
import picard.cmdline.argumentcollections.RequiredOutputArgumentCollection;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;


public class CollectMapqMetricsTest extends IntermediateFilesTest {
	public static Histogram<Integer> loadHistogram(File file) throws FileNotFoundException, IOException {
		try (FileReader reader  = new FileReader(file)) {
			MetricsFile<MapqMetrics,Integer> metricsFile = new MetricsFile<MapqMetrics, Integer>();
			metricsFile.read(reader);
			return metricsFile.getHistogram();
		}
	}
	@Test
	public void should_allow_no_histogram() throws IOException {
		CollectMapqMetrics cmm = new CollectMapqMetrics();
		cmm.INPUT = new File("src/test/resources/203541.bam");
		cmm.output = new RequiredOutputArgumentCollection(new File(testFolder.getRoot(), "mapqmetrics.txt"));
		SinglePassSamProgram.makeItSo(cmm.INPUT, null, true, 0, ImmutableList.of(cmm));
		assertTrue(cmm.output.getOutputFile().isFile());
	}
	@Test
	public void should_generate_metrics() throws IOException {
		CollectMapqMetrics cmm = new CollectMapqMetrics();
		cmm.INPUT = new File("src/test/resources/203541.bam");
		cmm.output = new RequiredOutputArgumentCollection(new File(testFolder.getRoot(), "mapqmetrics.txt"));
		cmm.Histogram_FILE = new File(testFolder.getRoot(), "mapqhistogram.pdf");
		SinglePassSamProgram.makeItSo(cmm.INPUT, null, true, 0, ImmutableList.of(cmm));
		assertTrue(cmm.output.getOutputFile().isFile());
	}
	@Test
	@Ignore("R dependency")
	public void should_generate_histogram() throws IOException {
		CollectMapqMetrics cmm = new CollectMapqMetrics();
		File histogram = new File(testFolder.getRoot(), "mapqhistogram.pdf");
		cmm.INPUT = new File("src/test/resources/203541.bam");
		cmm.output = new RequiredOutputArgumentCollection(new File(testFolder.getRoot(), "mapqmetrics.txt"));
		cmm.Histogram_FILE = histogram;
		SinglePassSamProgram.makeItSo(cmm.INPUT, null, true, 0, ImmutableList.of(cmm));
		assertTrue(histogram.isFile());
	}
	@Test
	public void metrics_should_correspond_to_read_counts() throws IOException {
		List<SAMRecord> reads = new ArrayList<>();
		for (int i = 0; i < 50; i++) {
			for (int j = 0; j <= i; j++) {
				reads.add(withMapq(j, Read(0, 1 + j,  10))[0]);
			}
		}
		createInput(reads);
		CollectMapqMetrics cmm = new CollectMapqMetrics();
		cmm.INPUT = input;
		cmm.output = new RequiredOutputArgumentCollection(new File(testFolder.getRoot(), "mapqmetrics.txt"));
		cmm.Histogram_FILE = new File(testFolder.getRoot(), "mapqhistogram.pdf");
		SinglePassSamProgram.makeItSo(cmm.INPUT, null, true, 0, ImmutableList.of(cmm));
		MapqMetrics metrics = IdsvSamFileMetrics.getMapqMetrics(cmm.output.getOutputFile());
		assertEquals(0, metrics.MIN_MAPQ);
		assertEquals(49, metrics.MAX_MAPQ);
		Histogram<Integer> histo = loadHistogram(cmm.output.getOutputFile());
		assertEquals(50, histo.values().size());
		for (Bin<Integer> bin : histo.values()) {
			assertEquals((int)bin.getId(), 50 - (int)bin.getValue());
		}
	}
}
