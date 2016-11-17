package au.edu.wehi.idsv.metrics;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.Test;

import com.google.common.collect.Lists;

import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import gridss.analysis.CigarDetailMetrics;
import gridss.analysis.IdsvMetrics;
import gridss.analysis.MapqMetrics;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import picard.analysis.InsertSizeMetrics;


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
		MetricsFile<CigarDetailMetrics, Integer> sc = new MetricsFile<CigarDetailMetrics, Integer>();
		MetricsFile<MapqMetrics, Integer> mq = new MetricsFile<MapqMetrics, Integer>();
		c.finish(is, idsv, mq, sc);
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
		MetricsFile<CigarDetailMetrics, Integer> sc = new MetricsFile<CigarDetailMetrics, Integer>();
		MetricsFile<MapqMetrics, Integer> mq = new MetricsFile<MapqMetrics, Integer>();
		c.finish(is, idsv, mq, sc);
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
		MetricsFile<CigarDetailMetrics, Integer> sc = new MetricsFile<CigarDetailMetrics, Integer>();
		MetricsFile<MapqMetrics, Integer> mq = new MetricsFile<MapqMetrics, Integer>();
		c.finish(is, idsv, mq, sc);
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
		MetricsFile<CigarDetailMetrics, Integer> sc = new MetricsFile<CigarDetailMetrics, Integer>();
		MetricsFile<MapqMetrics, Integer> mq = new MetricsFile<MapqMetrics, Integer>();
		c.finish(is, idsv, mq, sc);
		List<CigarDetailMetrics> msc = Lists.newArrayList(sc.getMetrics());
		msc.removeIf(cdm -> cdm.OPERATOR != CigarOperator.enumToCharacter(CigarOperator.SOFT_CLIP));
		assertEquals(6, msc.size());
		for (int i = 0; i < msc.size(); i++) {
			assertEquals(i, msc.get(i).LENGTH);
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
		MetricsFile<CigarDetailMetrics, Integer> sc = new MetricsFile<CigarDetailMetrics, Integer>();
		MetricsFile<MapqMetrics, Integer> mq = new MetricsFile<MapqMetrics, Integer>();
		c.finish(is, idsv, mq, sc);
		List<CigarDetailMetrics> msc = Lists.newArrayList(sc.getMetrics());
		msc.removeIf(cdm -> cdm.OPERATOR != CigarOperator.enumToCharacter(CigarOperator.SOFT_CLIP) || cdm.LENGTH < 1);
		assertEquals(3, msc.get(0).COUNT);
		assertEquals(0, msc.get(1).COUNT);
		assertEquals(0, msc.get(2).COUNT);
		assertEquals(0, msc.get(3).COUNT);
		assertEquals(1, msc.get(4).COUNT);
	}
	@Test
	public void should_add_implicit_zero_length_clips() {
		IdsvSamFileMetricsCollector c = new IdsvSamFileMetricsCollector(null);
		c.acceptRecord(Read(0, 1, "100M"), null); // 0, 0
		c.acceptRecord(Read(0, 1, "99M1S"), null); // 0, 1
		c.acceptRecord(Read(0, 1, "1S98M1S"), null);  // 1, 1
		c.acceptRecord(Read(0, 1, "95M5S"), null); // 0, 5
		MetricsFile<IdsvMetrics, Integer> idsv = new MetricsFile<IdsvMetrics, Integer>();
		MetricsFile<InsertSizeMetrics, Integer> is = new MetricsFile<InsertSizeMetrics, Integer>();
		MetricsFile<CigarDetailMetrics, Integer> sc = new MetricsFile<CigarDetailMetrics, Integer>();
		MetricsFile<MapqMetrics, Integer> mq = new MetricsFile<MapqMetrics, Integer>();
		c.finish(is, idsv, mq, sc);
		List<CigarDetailMetrics> msc = Lists.newArrayList(sc.getMetrics());
		msc.removeIf(cdm -> cdm.OPERATOR != CigarOperator.enumToCharacter(CigarOperator.SOFT_CLIP));
		assertEquals(4, msc.get(0).COUNT);
		msc = Lists.newArrayList(sc.getMetrics());
		msc.removeIf(cdm -> cdm.OPERATOR != CigarOperator.enumToCharacter(CigarOperator.HARD_CLIP));
		assertEquals(8, msc.get(0).COUNT);
	}
	@Test
	public void should_write_implicit_zero_cigar_for_missing_operators() {
		IdsvSamFileMetricsCollector c = new IdsvSamFileMetricsCollector(null);
		c.acceptRecord(Read(0, 1, "100M"), null); // 0, 0
		c.acceptRecord(Read(0, 1, "99M1S"), null); // 0, 1
		c.acceptRecord(Read(0, 1, "1S98M1S"), null);  // 1, 1
		c.acceptRecord(Read(0, 1, "95M5S"), null); // 0, 5
		MetricsFile<IdsvMetrics, Integer> idsv = new MetricsFile<IdsvMetrics, Integer>();
		MetricsFile<InsertSizeMetrics, Integer> is = new MetricsFile<InsertSizeMetrics, Integer>();
		MetricsFile<CigarDetailMetrics, Integer> sc = new MetricsFile<CigarDetailMetrics, Integer>();
		MetricsFile<MapqMetrics, Integer> mq = new MetricsFile<MapqMetrics, Integer>();
		c.finish(is, idsv, mq, sc);
		List<CigarDetailMetrics> msc = Lists.newArrayList(sc.getMetrics());
		msc.removeIf(cdm -> cdm.OPERATOR != CigarOperator.enumToCharacter(CigarOperator.X));
		assertEquals(4, msc.get(0).COUNT);
	}
	@Test
	public void should_calc_max_read_mapped_length() {
		IdsvSamFileMetricsCollector c = new IdsvSamFileMetricsCollector(null);
		c.acceptRecord(Read(0, 1, "5M90D5M100S"), null);
		MetricsFile<IdsvMetrics, Integer> idsv = new MetricsFile<IdsvMetrics, Integer>();
		MetricsFile<InsertSizeMetrics, Integer> is = new MetricsFile<InsertSizeMetrics, Integer>();
		MetricsFile<CigarDetailMetrics, Integer> sc = new MetricsFile<CigarDetailMetrics, Integer>();
		MetricsFile<MapqMetrics, Integer> mq = new MetricsFile<MapqMetrics, Integer>();
		c.finish(is, idsv, mq, sc);
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
		MetricsFile<CigarDetailMetrics, Integer> sc = new MetricsFile<CigarDetailMetrics, Integer>();
		MetricsFile<MapqMetrics, Integer> mq = new MetricsFile<MapqMetrics, Integer>();
		c.finish(is, idsv, mq, sc);
		assertEquals(11, (int)((IdsvMetrics)idsv.getMetrics().get(0)).MAX_PROPER_PAIR_FRAGMENT_LENGTH);
		
		c = new IdsvSamFileMetricsCollector(null);
		r = RP(0, 1, 100, 5)[1]; // mate before
		c.acceptRecord(r, null);
		c.finish(is, idsv, mq, sc);
		assertEquals(11, (int)((IdsvMetrics)idsv.getMetrics().get(0)).MAX_PROPER_PAIR_FRAGMENT_LENGTH);
	}
	@Test
	public void should_add_metrics_to_file() {
		IdsvSamFileMetricsCollector c = new IdsvSamFileMetricsCollector(null);
		MetricsFile<IdsvMetrics, Integer> idsv = new MetricsFile<IdsvMetrics, Integer>();
		MetricsFile<InsertSizeMetrics, Integer> is = new MetricsFile<InsertSizeMetrics, Integer>();
		MetricsFile<CigarDetailMetrics, Integer> sc = new MetricsFile<CigarDetailMetrics, Integer>();
		MetricsFile<MapqMetrics, Integer> mq = new MetricsFile<MapqMetrics, Integer>();
		c.finish(is, idsv, mq, sc);
		assertEquals(1, idsv.getMetrics().size());
		assertTrue(idsv.getMetrics().get(0) instanceof IdsvMetrics);
	}
	@Test
	public void fragment_size_should_not_be_set_if_no_proper_pairs_found() {
		IdsvSamFileMetricsCollector c = new IdsvSamFileMetricsCollector(null);
		c.acceptRecord(Read(0, 1, "100M"), null);
		MetricsFile<IdsvMetrics, Integer> idsv = new MetricsFile<IdsvMetrics, Integer>();
		MetricsFile<InsertSizeMetrics, Integer> is = new MetricsFile<InsertSizeMetrics, Integer>();
		MetricsFile<CigarDetailMetrics, Integer> sc = new MetricsFile<CigarDetailMetrics, Integer>();
		MetricsFile<MapqMetrics, Integer> mq = new MetricsFile<MapqMetrics, Integer>();
		c.finish(is, idsv, mq, sc);
		assertNull(idsv.getMetrics().get(0).MAX_PROPER_PAIR_FRAGMENT_LENGTH);
		assertNull(idsv.getMetrics().get(0).MIN_PROPER_PAIR_FRAGMENT_LENGTH);
	}
	@Test
	public void mapq_should_be_calculated() {
		IdsvSamFileMetricsCollector c = new IdsvSamFileMetricsCollector(null);
		c.acceptRecord(withMapq(2, Read(0, 1, "100M"))[0], null);
		MetricsFile<IdsvMetrics, Integer> idsv = new MetricsFile<IdsvMetrics, Integer>();
		MetricsFile<InsertSizeMetrics, Integer> is = new MetricsFile<InsertSizeMetrics, Integer>();
		MetricsFile<CigarDetailMetrics, Integer> sc = new MetricsFile<CigarDetailMetrics, Integer>();
		MetricsFile<MapqMetrics, Integer> mq = new MetricsFile<MapqMetrics, Integer>();
		c.finish(is, idsv, mq, sc);
		assertEquals(2, mq.getMetrics().get(0).MIN_MAPQ);
		assertEquals(2, mq.getMetrics().get(0).MAX_MAPQ);
	}
	@Test
	public void mapq_should_be_extrema() {
		IdsvSamFileMetricsCollector c = new IdsvSamFileMetricsCollector(null);
		c.acceptRecord(withMapq(1, Read(0, 1, "100M"))[0], null);
		c.acceptRecord(withMapq(2, Read(0, 1, "101M"))[0], null);
		c.acceptRecord(withMapq(2, Read(0, 1, "102M"))[0], null);
		c.acceptRecord(withMapq(3, Read(0, 1, "103M"))[0], null);
		MetricsFile<IdsvMetrics, Integer> idsv = new MetricsFile<IdsvMetrics, Integer>();
		MetricsFile<InsertSizeMetrics, Integer> is = new MetricsFile<InsertSizeMetrics, Integer>();
		MetricsFile<CigarDetailMetrics, Integer> sc = new MetricsFile<CigarDetailMetrics, Integer>();
		MetricsFile<MapqMetrics, Integer> mq = new MetricsFile<MapqMetrics, Integer>();
		c.finish(is, idsv, mq, sc);
		assertEquals(1, mq.getMetrics().get(0).MIN_MAPQ);
		assertEquals(3, mq.getMetrics().get(0).MAX_MAPQ);
	}
}
