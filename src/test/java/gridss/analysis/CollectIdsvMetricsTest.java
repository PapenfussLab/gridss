package gridss.analysis;

import static org.junit.Assert.assertEquals;

import java.io.File;

import org.junit.Test;

import au.edu.wehi.idsv.IntermediateFilesTest;
import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.SAMRecord;

public class CollectIdsvMetricsTest extends IntermediateFilesTest {
	@Test
	public void should_calc_max_read_length() {
		File metricFiles = new File(testFolder.getRoot(), "metrics.txt");
		CollectIdsvMetrics c = new CollectIdsvMetrics();
		c.OUTPUT = metricFiles;
		c.setup(null, null);
		c.acceptRead(RP(0, 1, 2, 1)[0], null);
		c.acceptRead(RP(0, 1, 2, 1)[1], null);
		c.acceptRead(RP(0, 1, 7, 5)[0], null);
		c.acceptRead(RP(0, 1, 7, 5)[1], null);
		c.acceptRead(Read(0, 1, "100M"), null);
		c.finish();
		IdsvMetrics metrics = new IdsvSamFileMetrics(null, metricFiles, null, null, false).getIdsvMetrics(); 
		assertEquals(100, metrics.MAX_READ_LENGTH);
	}
	@Test
	public void should_calc_read_statistics() {
		File metricFiles = new File(testFolder.getRoot(), "metrics.txt");
		CollectIdsvMetrics c = new CollectIdsvMetrics();
		c.OUTPUT = metricFiles;
		c.setup(null, null);
		c.acceptRead(RP(0, 1, 2, 1)[1], null);
		c.acceptRead(RP(0, 1, 2, 1)[0], null);
		c.acceptRead(RP(0, 1, 7, 5)[0], null);
		c.acceptRead(RP(0, 1, 7, 5)[1], null);
		c.acceptRead(Read(0, 1, "100M"), null);
		c.acceptRead(OEA(0, 1, "1M", true)[0], null);
		c.acceptRead(OEA(0, 1, "1M", true)[1], null);
		c.acceptRead(Unmapped(100), null);
		SAMRecord[] unmapped = RP(0, 1, 2, 1);
		unmapped[0].setReadUnmappedFlag(true);
		unmapped[1].setReadUnmappedFlag(true);
		SAMRecordUtil.pairReads(unmapped[0], unmapped[1]);
		c.acceptRead(unmapped[0], null);
		c.acceptRead(unmapped[1], null);
		c.finish();
		IdsvMetrics metrics = new IdsvSamFileMetrics(null, metricFiles, null, null, false).getIdsvMetrics(); 
		assertEquals(10, metrics.READS);
		assertEquals(6, metrics.MAPPED_READS);
	}
	@Test
	public void should_calc_read_pairing_statistics() {
		File metricFiles = new File(testFolder.getRoot(), "metrics.txt");
		CollectIdsvMetrics c = new CollectIdsvMetrics();
		c.OUTPUT = metricFiles;
		c.setup(null, null);
		c.acceptRead(RP(0, 1, 2, 1)[1], null);
		c.acceptRead(RP(0, 1, 2, 1)[0], null);
		c.acceptRead(RP(0, 1, 7, 5)[0], null);
		c.acceptRead(RP(0, 1, 7, 5)[1], null);
		c.acceptRead(Read(0, 1, "100M"), null);
		c.acceptRead(OEA(0, 1, "1M", true)[0], null);
		c.acceptRead(OEA(0, 1, "1M", true)[1], null);
		c.acceptRead(Unmapped(100), null);
		SAMRecord[] unmapped = RP(0, 1, 2, 1);
		unmapped[0].setReadUnmappedFlag(true);
		unmapped[1].setReadUnmappedFlag(true);
		SAMRecordUtil.pairReads(unmapped[0], unmapped[1]);
		c.acceptRead(unmapped[0], null);
		c.acceptRead(unmapped[1], null);
		c.finish();
		IdsvMetrics metrics = new IdsvSamFileMetrics(null, metricFiles, null, null, false).getIdsvMetrics();
		assertEquals(4, metrics.READ_PAIRS);
		assertEquals(2, metrics.READ_PAIRS_BOTH_MAPPED);
		assertEquals(1, metrics.READ_PAIRS_ONE_MAPPED);
		assertEquals(1, metrics.READ_PAIRS_ZERO_MAPPED);
	}
	@Test
	public void should_calc_alternate_mappings() {
		File metricFiles = new File(testFolder.getRoot(), "metrics.txt");
		CollectIdsvMetrics c = new CollectIdsvMetrics();
		c.OUTPUT = metricFiles;
		c.setup(null, null);
		c.acceptRead(Read(0, 1, "1M"), null);
		SAMRecord r2 = Read(0, 1, "1M");
		r2.setSupplementaryAlignmentFlag(true);
		c.acceptRead(r2, null);
		SAMRecord r3 = Read(0, 1, "1M");
		r3.setSecondaryAlignment(true);
		c.acceptRead(r3, null);
		SAMRecord r4 = Read(0, 1, "1M");
		r4.setSecondaryAlignment(true);
		r4.setAttribute("SA", "chr1,5,+,1M,0,0");
		c.acceptRead(r4, null);
		c.finish();
		IdsvMetrics metrics = new IdsvSamFileMetrics(null, metricFiles, null, null, false).getIdsvMetrics();
		assertEquals(1, metrics.SECONDARY_NOT_SPLIT);
	}
	@Test
	public void should_not_count_secondary() {
		File metricFiles = new File(testFolder.getRoot(), "metrics.txt");
		CollectIdsvMetrics c = new CollectIdsvMetrics();
		c.OUTPUT = metricFiles;
		c.setup(null, null);
		SAMRecord r1 = Read(0, 1, "1M");
		r1.setSecondaryAlignment(true);
		c.acceptRead(r1, null);
		c.finish();
		IdsvMetrics metrics = new IdsvSamFileMetrics(null, metricFiles, null, null, false).getIdsvMetrics();
		assertEquals(0, metrics.READS);
		
		c = new CollectIdsvMetrics();
		c.OUTPUT = metricFiles;
		c.COUNT_SECONDARY = true;
		c.setup(null, null);
		c.acceptRead(r1, null);
		c.finish();
		IdsvMetrics metrics_counting = new IdsvSamFileMetrics(null, metricFiles, null, null, false).getIdsvMetrics();
		assertEquals(1, metrics_counting.READS);
	}
	@Test
	public void should_not_count_supplementary() {
		File metricFiles = new File(testFolder.getRoot(), "metrics.txt");
		CollectIdsvMetrics c = new CollectIdsvMetrics();
		c.OUTPUT = metricFiles;
		c.setup(null, null);
		SAMRecord r1 = Read(0, 1, "1M");
		r1.setSupplementaryAlignmentFlag(true);
		c.acceptRead(r1, null);
		c.finish();
		IdsvMetrics metrics = new IdsvSamFileMetrics(null, metricFiles, null, null, false).getIdsvMetrics();
		assertEquals(0, metrics.READS);
	}
}
