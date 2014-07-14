package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;

import java.io.File;
import java.io.IOException;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import picard.analysis.InsertSizeMetrics;
import picard.analysis.directed.InsertSizeMetricsCollector;

public class PicardMetricsTest extends TestHelper {
	@Rule
    public TemporaryFolder testFolder = new TemporaryFolder();
	private InsertSizeMetricsCollector readPair(int start1, int start2, int readLength) {
		InsertSizeMetricsCollector x = PicardMetrics.createCollector(null);
		// create read pair
		SAMRecord[] pair = RP(0, start1, start2, readLength);
		x.acceptRecord(pair[0], null);
		x.acceptRecord(pair[1], null);
		x.finish();
		return x;
	}
	private InsertSizeMetricsCollector readCollector(SAMRecord[]... args) {
		InsertSizeMetricsCollector x = PicardMetrics.createCollector(null);
		for (SAMRecord[] ra : args) {
			for (SAMRecord r : ra) {
				x.acceptRecord(r, null);
			}
		}
		x.finish();
		return x;
	}
	@Test
	public void shouldSaveMetricsFile() throws IOException {
		InsertSizeMetricsCollector x = readPair(1, 2, 1);
		File f = testFolder.newFile();
		f.delete();
		PicardMetrics.save(x, new MetricsFile<InsertSizeMetrics, Integer>(), f);
		f = new File(f.toString());
		assertTrue(f.exists());
	}
	@Test
	public void shouldLoadMetricsFromFile() throws IOException {
		int readLength = 50;
		File f = testFolder.newFile();
		PicardMetrics.save(readPair(100, 200, readLength), new MetricsFile<InsertSizeMetrics, Integer>(), f);
		RelevantMetrics metrics = new PicardMetrics(f);
		assertEquals(150, metrics.getMedianFragmentSize(), 0);
	}
	@Test
	public void shouldLoadMetricsFromFile_with_no_reads() throws IOException {
		File f = testFolder.newFile();
		InsertSizeMetricsCollector x = PicardMetrics.createCollector(null);
		x.finish();
		PicardMetrics.save(x, new MetricsFile<InsertSizeMetrics, Integer>(), f);
		RelevantMetrics metrics = new PicardMetrics(f);
		assertEquals(0, metrics.getMedianFragmentSize(), 0);
		assertEquals(0, metrics.getFragmentSizeStdDev(), 0);
		assertEquals(0, metrics.getMaxFragmentSize());
	}
	@Test
	public void getMedianFragmentSize_should_return_median_fragment_size() throws IOException {
		int readLength = 100;
		File f = testFolder.newFile();
		PicardMetrics.save(readPair(100, 200, readLength), new MetricsFile<InsertSizeMetrics, Integer>(), f);
		RelevantMetrics metrics = new PicardMetrics(f);
		assertEquals(100 + readLength, metrics.getMedianFragmentSize(), 0);
		
		readLength = 50;
		 f = testFolder.newFile();
		PicardMetrics.save(readPair(100, 200, readLength), new MetricsFile<InsertSizeMetrics, Integer>(), f);
		metrics = new PicardMetrics(f);
		assertEquals(100 + readLength, metrics.getMedianFragmentSize(), 0);
	}
	@Test
	public void shouldUseMADforStdDev() throws IOException {
		// 3 pairs
		// (100,300) (100, 250), (100, 350) 
		File f = testFolder.newFile();
		PicardMetrics.save(readCollector(
				RP(0, 100, 300, 100),
				RP(0, 100, 250, 100),
				RP(0, 100, 350, 100)
				), new MetricsFile<InsertSizeMetrics, Integer>(), f);
		RelevantMetrics metrics = new PicardMetrics(f);
		assertEquals(1.4826 * 50, metrics.getFragmentSizeStdDev(), 0.001);
	}
}
