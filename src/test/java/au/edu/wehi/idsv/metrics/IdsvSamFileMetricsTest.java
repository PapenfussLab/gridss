package au.edu.wehi.idsv.metrics;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import picard.analysis.InsertSizeMetrics;
import au.edu.wehi.idsv.TestHelper;

public class IdsvSamFileMetricsTest extends TestHelper {
	@Test
	public void wrapper_inner_metrics() {
		IdsvMetrics im = new IdsvMetrics();
		InsertSizeMetrics ism = new InsertSizeMetrics();
		List<SoftClipDetailMetrics> sc = new ArrayList<SoftClipDetailMetrics>();
		InsertSizeDistribution isd = new InsertSizeDistribution(new int[] { 1}, new double[] { 1}, 1);
		IdsvSamFileMetrics metrics = new IdsvSamFileMetrics(ism, im, isd, sc);
		assertEquals(im, metrics.getIdsvMetrics());
		assertEquals(isd, metrics.getInsertSizeDistribution());
		assertEquals(ism, metrics.getInsertSizeMetrics());
		assertEquals(sc, metrics.getSoftClipDetailMetrics());
	}
	/*
	@Test
	public void shouldUseMADforStdDev() {
		IdsvMetrics im = new IdsvMetrics();
		InsertSizeMetrics ism = new InsertSizeMetrics();
		ism.MEDIAN_ABSOLUTE_DEVIATION = 10;
		IdsvSamFileMetrics metrics = new IdsvSamFileMetrics(getContext(), ism, im, null);
		assertEquals(10 * 1.4826, metrics.getFragmentSizeStdDev(), 0.0001);
	}
	@Test
	public void shouldUseMaxProperPairSizeForMaxFragmentSize() {
		IdsvMetrics im = new IdsvMetrics();
		im.MAX_PROPER_PAIR_FRAGMENT_LENGTH = 10;
		InsertSizeMetrics ism = new InsertSizeMetrics();
		IdsvSamFileMetrics metrics = new IdsvSamFileMetrics(getContext(), ism, im, null);
		assertEquals(10, metrics.getMaxFragmentSize());
	}
	@Test
	public void should_use_read_length_as_fragment_size_for_unpaired_reads() {
		IdsvMetrics im = new IdsvMetrics();
		im.MAX_PROPER_PAIR_FRAGMENT_LENGTH = 10;
		im.MAX_READ_LENGTH = 50;
		InsertSizeMetrics ism = new InsertSizeMetrics();
		IdsvSamFileMetrics metrics = new IdsvSamFileMetrics(getContext(), ism, im, null);
		assertEquals(50, metrics.getMaxFragmentSize());
	}
	@Test
	public void should_use_concordant_percentile_for_max_fragment_size_if_not_using_proper_pair_flag() {
		ProcessingContext pc = getContext();
		pc.getReadPairParameters().useProperPairFlag = false;
		pc.getReadPairParameters().concordantPercent = 0.5;
		IdsvMetrics im = new IdsvMetrics();
		im.MAX_PROPER_PAIR_FRAGMENT_LENGTH = 10;
		im.MAX_READ_LENGTH = 50;
		InsertSizeMetrics ism = new InsertSizeMetrics();
		int[] fragSize = new int[1000];
		double[] probability = new double[1000];
		for (int i = 0; i < 1000; i++) {
			fragSize[i] = i;
			probability[i] = 1d / 1000d;
		}
		InsertSizeDistribution isd = new InsertSizeDistribution(fragSize, probability, 10000);
		IdsvSamFileMetrics metrics = new IdsvSamFileMetrics(pc, ism, im, isd);
		assertEquals(749, metrics.getMaxFragmentSize());
	}
	*/
}
