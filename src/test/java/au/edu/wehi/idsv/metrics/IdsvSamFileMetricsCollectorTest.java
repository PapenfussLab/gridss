package au.edu.wehi.idsv.metrics;

import static org.junit.Assert.*;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;

import org.junit.Test;

import picard.analysis.InsertSizeMetrics;
import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.sam.SAMRecordUtil;


public class IdsvSamFileMetricsCollectorTest extends TestHelper {
	@Test
	public void should_calc_max_read_length() {
		IdsvSamFileMetricsCollector c = new IdsvSamFileMetricsCollector(null);
		c.acceptRecord(RP(0, 1, 2, 1)[0], null);
		c.acceptRecord(RP(0, 1, 2, 1)[1], null);
		c.acceptRecord(RP(0, 1, 7, 5)[0], null);
		c.acceptRecord(RP(0, 1, 7, 5)[1], null);
		c.acceptRecord(Read(0, 1, "100M"), null);
		MetricsFile<IdsvMetrics, Integer> idsv = new MetricsFile<IdsvMetrics, Integer>();
		MetricsFile<InsertSizeMetrics, Integer> is = new MetricsFile<InsertSizeMetrics, Integer>();
		MetricsFile<SoftClipDetailMetrics, Integer> sc = new MetricsFile<SoftClipDetailMetrics, Integer>();
		c.finish(is, idsv, sc);
		assertEquals(100, (int)((IdsvMetrics)idsv.getMetrics().get(0)).MAX_READ_LENGTH);
	}
	@Test
	public void should_calc_read_statistics() {
		IdsvSamFileMetricsCollector c = new IdsvSamFileMetricsCollector(null);
		c.acceptRecord(RP(0, 1, 2, 1)[1], null);
		c.acceptRecord(RP(0, 1, 2, 1)[0], null);
		c.acceptRecord(RP(0, 1, 7, 5)[0], null);
		c.acceptRecord(RP(0, 1, 7, 5)[1], null);
		c.acceptRecord(Read(0, 1, "100M"), null);
		c.acceptRecord(OEA(0, 1, "1M", true)[0], null);
		c.acceptRecord(OEA(0, 1, "1M", true)[1], null);
		c.acceptRecord(Unmapped(100), null);
		SAMRecord[] unmapped = RP(0, 1, 2, 1);
		unmapped[0].setReadUnmappedFlag(true);
		unmapped[1].setReadUnmappedFlag(true);
		SAMRecordUtil.pairReads(unmapped[0], unmapped[1]);
		c.acceptRecord(unmapped[0], null);
		c.acceptRecord(unmapped[1], null);
		MetricsFile<IdsvMetrics, Integer> idsv = new MetricsFile<IdsvMetrics, Integer>();
		MetricsFile<InsertSizeMetrics, Integer> is = new MetricsFile<InsertSizeMetrics, Integer>();
		MetricsFile<SoftClipDetailMetrics, Integer> sc = new MetricsFile<SoftClipDetailMetrics, Integer>();
		c.finish(is, idsv, sc);
		assertEquals(10, idsv.getMetrics().get(0).READS);
		assertEquals(6, idsv.getMetrics().get(0).MAPPED_READS);
	}
	@Test
	public void should_calc_read_pairing_statistics() {
		IdsvSamFileMetricsCollector c = new IdsvSamFileMetricsCollector(null);
		c.acceptRecord(RP(0, 1, 2, 1)[1], null);
		c.acceptRecord(RP(0, 1, 2, 1)[0], null);
		c.acceptRecord(RP(0, 1, 7, 5)[0], null);
		c.acceptRecord(RP(0, 1, 7, 5)[1], null);
		c.acceptRecord(Read(0, 1, "100M"), null);
		c.acceptRecord(OEA(0, 1, "1M", true)[0], null);
		c.acceptRecord(OEA(0, 1, "1M", true)[1], null);
		c.acceptRecord(Unmapped(100), null);
		SAMRecord[] unmapped = RP(0, 1, 2, 1);
		unmapped[0].setReadUnmappedFlag(true);
		unmapped[1].setReadUnmappedFlag(true);
		SAMRecordUtil.pairReads(unmapped[0], unmapped[1]);
		c.acceptRecord(unmapped[0], null);
		c.acceptRecord(unmapped[1], null);
		MetricsFile<IdsvMetrics, Integer> idsv = new MetricsFile<IdsvMetrics, Integer>();
		MetricsFile<InsertSizeMetrics, Integer> is = new MetricsFile<InsertSizeMetrics, Integer>();
		MetricsFile<SoftClipDetailMetrics, Integer> sc = new MetricsFile<SoftClipDetailMetrics, Integer>();
		c.finish(is, idsv, sc);
		assertEquals(4, idsv.getMetrics().get(0).READ_PAIRS);
		assertEquals(2, idsv.getMetrics().get(0).READ_PAIRS_BOTH_MAPPED);
		assertEquals(1, idsv.getMetrics().get(0).READ_PAIRS_ONE_MAPPED);
		assertEquals(1, idsv.getMetrics().get(0).READ_PAIRS_ZERO_MAPPED);
	}
	@Test
	public void should_create_soft_clip_metrics_up_to_read_length() {
		IdsvSamFileMetricsCollector c = new IdsvSamFileMetricsCollector(null);
		c.acceptRecord(Read(0, 1, "100M"), null);
		c.acceptRecord(Read(0, 1, "99M1S"), null);
		c.acceptRecord(Read(0, 1, "1S98M1S"), null);
		c.acceptRecord(Read(0, 1, "95M5S"), null);
		MetricsFile<IdsvMetrics, Integer> idsv = new MetricsFile<IdsvMetrics, Integer>();
		MetricsFile<InsertSizeMetrics, Integer> is = new MetricsFile<InsertSizeMetrics, Integer>();
		MetricsFile<SoftClipDetailMetrics, Integer> sc = new MetricsFile<SoftClipDetailMetrics, Integer>();
		c.finish(is, idsv, sc);
		assertEquals(100, sc.getMetrics().size());
		for (int i = 0; i < sc.getMetrics().size(); i++) {
			assertEquals(i, sc.getMetrics().get(i).LENGTH);
		}
	}
	@Test
	public void should_count_soft_clips_on_both_ends() {
		IdsvSamFileMetricsCollector c = new IdsvSamFileMetricsCollector(null);
		c.acceptRecord(Read(0, 1, "100M"), null); // 0, 0
		c.acceptRecord(Read(0, 1, "99M1S"), null); // 0, 1
		c.acceptRecord(Read(0, 1, "1S98M1S"), null);  // 1, 1
		c.acceptRecord(Read(0, 1, "95M5S"), null); // 0, 5
		MetricsFile<IdsvMetrics, Integer> idsv = new MetricsFile<IdsvMetrics, Integer>();
		MetricsFile<InsertSizeMetrics, Integer> is = new MetricsFile<InsertSizeMetrics, Integer>();
		MetricsFile<SoftClipDetailMetrics, Integer> sc = new MetricsFile<SoftClipDetailMetrics, Integer>();
		c.finish(is, idsv, sc);
		assertEquals(4, ((SoftClipDetailMetrics)sc.getMetrics().get(0)).READCOUNT);
		assertEquals(3, ((SoftClipDetailMetrics)sc.getMetrics().get(1)).READCOUNT);
		assertEquals(0, ((SoftClipDetailMetrics)sc.getMetrics().get(2)).READCOUNT);
		assertEquals(0, ((SoftClipDetailMetrics)sc.getMetrics().get(3)).READCOUNT);
		assertEquals(0, ((SoftClipDetailMetrics)sc.getMetrics().get(4)).READCOUNT);
		assertEquals(1, ((SoftClipDetailMetrics)sc.getMetrics().get(5)).READCOUNT);
		assertEquals(0, ((SoftClipDetailMetrics)sc.getMetrics().get(6)).READCOUNT);
	}
	@Test
	public void should_calc_max_read_mapped_length() {
		IdsvSamFileMetricsCollector c = new IdsvSamFileMetricsCollector(null);
		c.acceptRecord(Read(0, 1, "5M90D5M100S"), null);
		MetricsFile<IdsvMetrics, Integer> idsv = new MetricsFile<IdsvMetrics, Integer>();
		MetricsFile<InsertSizeMetrics, Integer> is = new MetricsFile<InsertSizeMetrics, Integer>();
		MetricsFile<SoftClipDetailMetrics, Integer> sc = new MetricsFile<SoftClipDetailMetrics, Integer>();
		c.finish(is, idsv, sc);
		assertEquals(100, (int)((IdsvMetrics)idsv.getMetrics().get(0)).MAX_READ_MAPPED_LENGTH);
	}
	@Test
	public void should_calc_max_fragment_size() {
		IdsvSamFileMetricsCollector c = new IdsvSamFileMetricsCollector(null);
		c.acceptRecord(Read(0, 1, "100M"), null);
		c.acceptRecord(RP(0, 1, 2, 1)[0], null);
		c.acceptRecord(RP(0, 1, 7, 5)[0], null);
		SAMRecord r = RP(0, 1, 100, 5)[0]; // should ignore this as it's not a proper pair
		r.setProperPairFlag(false);
		c.acceptRecord(r, null);
		
		// 12345678901234567890
		// ----> <----
		c.acceptRecord(Read(0, 1, "100M"), null);
		MetricsFile<IdsvMetrics, Integer> idsv = new MetricsFile<IdsvMetrics, Integer>();
		MetricsFile<InsertSizeMetrics, Integer> is = new MetricsFile<InsertSizeMetrics, Integer>();
		MetricsFile<SoftClipDetailMetrics, Integer> sc = new MetricsFile<SoftClipDetailMetrics, Integer>();
		c.finish(is, idsv, sc);
		assertEquals(11, (int)((IdsvMetrics)idsv.getMetrics().get(0)).MAX_PROPER_PAIR_FRAGMENT_LENGTH);
		
		c = new IdsvSamFileMetricsCollector(null);
		r = RP(0, 1, 100, 5)[1]; // mate before
		c.acceptRecord(r, null);
		c.finish(is, idsv, sc);
		assertEquals(11, (int)((IdsvMetrics)idsv.getMetrics().get(0)).MAX_PROPER_PAIR_FRAGMENT_LENGTH);
	}
	@Test
	public void should_add_metrics_to_file() {
		IdsvSamFileMetricsCollector c = new IdsvSamFileMetricsCollector(null);
		MetricsFile<IdsvMetrics, Integer> idsv = new MetricsFile<IdsvMetrics, Integer>();
		MetricsFile<InsertSizeMetrics, Integer> is = new MetricsFile<InsertSizeMetrics, Integer>();
		MetricsFile<SoftClipDetailMetrics, Integer> sc = new MetricsFile<SoftClipDetailMetrics, Integer>();
		c.finish(is, idsv, sc);
		assertEquals(1, idsv.getMetrics().size());
		assertTrue(idsv.getMetrics().get(0) instanceof IdsvMetrics);
	}
	@Test
	public void fragment_size_should_not_be_set_if_no_proper_pairs_found() {
		IdsvSamFileMetricsCollector c = new IdsvSamFileMetricsCollector(null);
		c.acceptRecord(Read(0, 1, "100M"), null);
		MetricsFile<IdsvMetrics, Integer> idsv = new MetricsFile<IdsvMetrics, Integer>();
		MetricsFile<InsertSizeMetrics, Integer> is = new MetricsFile<InsertSizeMetrics, Integer>();
		MetricsFile<SoftClipDetailMetrics, Integer> sc = new MetricsFile<SoftClipDetailMetrics, Integer>();
		c.finish(is, idsv, sc);
		assertNull(idsv.getMetrics().get(0).MAX_PROPER_PAIR_FRAGMENT_LENGTH);
		assertNull(idsv.getMetrics().get(0).MIN_PROPER_PAIR_FRAGMENT_LENGTH);
	}
}
