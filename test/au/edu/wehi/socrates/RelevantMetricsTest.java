package au.edu.wehi.socrates;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;

import net.sf.picard.analysis.InsertSizeMetrics;
import net.sf.picard.analysis.directed.InsertSizeMetricsCollector;
import net.sf.picard.metrics.MetricsFile;
import net.sf.samtools.SAMRecord;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

public class RelevantMetricsTest {
	@Rule
    public TemporaryFolder testFolder = new TemporaryFolder();
	private InsertSizeMetricsCollector readPair(int start1, int start2, int readLength) {
		InsertSizeMetricsCollector x = RelevantMetrics.createCollector(null);
		// create read pair
		SAMRecord r1, r2;
		x.acceptRecord(r1, null);
		x.acceptRecord(r2, null);
		return x;
	}
	@Test
	public void shouldSaveMetricsFile() throws IOException {
		InsertSizeMetricsCollector x = readPair(1, 2, 1);
		File f = testFolder.newFile();
		f.delete();
		RelevantMetrics.save(x, new MetricsFile<InsertSizeMetrics, Integer>(), f);
		f = new File(f.toString());
		assertTrue(f.exists());
	}
	@Test
	public void shouldLoadMetricsFromFile() throws IOException {
		File f = testFolder.newFile();
		RelevantMetrics.save(readPair(100, 200, 100), new MetricsFile<InsertSizeMetrics, Integer>(), f);
		RelevantMetrics metrics = new RelevantMetrics(f);
		assertEquals(300, metrics.getMedianFragmentSize(), 0);
	}
	@Test
	public void shouldUseMADforStdDev() throws IOException {
		// 3 pairs
		// (100,300) (100, 250), (100, 350) 
		File f = testFolder.newFile();
		RelevantMetrics.save(readPair(), new MetricsFile<InsertSizeMetrics, Integer>(), f);
		RelevantMetrics metrics = new RelevantMetrics(f);
		assertEquals(1.4826 * 50, metrics.getFragmentSizeStdDev(), 0.001);
	}
}
