package gridss.analysis;

import static org.junit.Assert.assertEquals;

import java.io.File;

import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;

import au.edu.wehi.idsv.IntermediateFilesTest;
import au.edu.wehi.idsv.ReadPairConcordanceMethod;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.metrics.MetricsFile;

public class CollectStructuralVariantReadMetricsTest extends IntermediateFilesTest {
	private static <T> T getMetrics(File file, Class<T> clazz) {
		for (T metric : Iterables.filter(MetricsFile.readBeans(file), clazz)) {
			return metric;
		}
		return null;
	}
	@Test
	public void should_consider_multiple_discordant_pairs() {
		File metricFiles = new File(testFolder.getRoot(), "metrics.txt");
		CollectStructuralVariantReadMetrics c = new CollectStructuralVariantReadMetrics();
		c.OUTPUT = metricFiles;
		c.READ_PAIR_CONCORDANCE_METHOD = ReadPairConcordanceMethod.FIXED;
		c.FIXED_READ_PAIR_CONCORDANCE_MAX_FRAGMENT_SIZE = 1000;
		c.FIXED_READ_PAIR_CONCORDANCE_MIN_FRAGMENT_SIZE = 1000;
		c.setup(null, null);
		c.acceptFragment(ImmutableList.of(
				withReadName("read1", RP(0, 1, 2, 1))[0],
				withReadName("read1", RP(1, 1, 2, 1))[0],
				withReadName("read1", RP(0, 1, 2, 1))[1],
				withReadName("read1", RP(1, 1, 2, 1))[1]), null);
		c.finish();
		StructuralVariantReadMetrics metrics = getMetrics(metricFiles, StructuralVariantReadMetrics.class);
		assertEquals(1, metrics.DISCORDANT_READ_PAIRS);
		assertEquals(4, metrics.DISCORDANT_READ_PAIR_ALIGNMENTS);
		assertEquals(1, metrics.STRUCTURAL_VARIANT_READ_PAIRS);
		assertEquals(0, metrics.UNMAPPED_MATE_READ_ALIGNMENTS);
		assertEquals(0, metrics.UNMAPPED_MATE_READS);
	}
	@Test
	public void should_consider_multiple_oeas() {
		File metricFiles = new File(testFolder.getRoot(), "metrics.txt");
		CollectStructuralVariantReadMetrics c = new CollectStructuralVariantReadMetrics();
		c.OUTPUT = metricFiles;
		c.READ_PAIR_CONCORDANCE_METHOD = ReadPairConcordanceMethod.FIXED;
		c.FIXED_READ_PAIR_CONCORDANCE_MAX_FRAGMENT_SIZE = 1000;
		c.FIXED_READ_PAIR_CONCORDANCE_MIN_FRAGMENT_SIZE = 1000;
		c.setup(null, null);
		c.acceptFragment(ImmutableList.of(
				withReadName("read1", OEA(0, 1, "1M", true))[0],
				withReadName("read1", OEA(0, 1, "1M", true))[1],
				withReadName("read1", OEA(0, 2, "1M", true))[0],
				withReadName("read1", OEA(0, 2, "1M", true))[1]), null);
		c.finish();
		StructuralVariantReadMetrics metrics = getMetrics(metricFiles, StructuralVariantReadMetrics.class);
		assertEquals(0, metrics.DISCORDANT_READ_PAIRS);
		assertEquals(0, metrics.DISCORDANT_READ_PAIR_ALIGNMENTS);
		assertEquals(1, metrics.STRUCTURAL_VARIANT_READ_PAIRS);
		assertEquals(2, metrics.UNMAPPED_MATE_READ_ALIGNMENTS);
		assertEquals(1, metrics.UNMAPPED_MATE_READS);
	}
	@Test
	public void should_count_indel_alignments_and_reads_separately() {
		File metricFiles = new File(testFolder.getRoot(), "metrics.txt");
		CollectStructuralVariantReadMetrics c = new CollectStructuralVariantReadMetrics();
		c.OUTPUT = metricFiles;
		c.setup(null, null);
		c.acceptFragment(ImmutableList.of(
				withReadName("indel", Read(0, 1, "10M5I10M"))[0],
				withReadName("indel", Read(0, 2, "10M5I10M"))[0]), null);
		c.finish();
		StructuralVariantReadMetrics metrics = getMetrics(metricFiles, StructuralVariantReadMetrics.class);
		assertEquals(1, metrics.INDEL_READS);
		assertEquals(2, metrics.INDEL_READ_ALIGNMENTS);
		assertEquals(1, metrics.STRUCTURAL_VARIANT_READS);
		assertEquals(2, metrics.STRUCTURAL_VARIANT_READ_ALIGNMENTS);
	}
	@Test
	public void should_count_sc_alignments_and_reads_separately() {
		File metricFiles = new File(testFolder.getRoot(), "metrics.txt");
		CollectStructuralVariantReadMetrics c = new CollectStructuralVariantReadMetrics();
		c.OUTPUT = metricFiles;
		c.setup(null, null);
		c.acceptFragment(ImmutableList.of(
				withReadName("sc", Read(0, 1, "10M5S"))[0],
				withReadName("sc", Read(0, 2, "5S10M5S"))[0]), null);
		c.finish();
		StructuralVariantReadMetrics metrics = getMetrics(metricFiles, StructuralVariantReadMetrics.class);
		assertEquals(1, metrics.SOFT_CLIPPED_READS);
		assertEquals(2, metrics.SOFT_CLIPPED_READ_ALIGNMENTS);
		assertEquals(1, metrics.STRUCTURAL_VARIANT_READS);
		assertEquals(2, metrics.STRUCTURAL_VARIANT_READ_ALIGNMENTS);
	}
	@Test
	public void should_iterate_by_query_name() {
		File metricFiles = new File(testFolder.getRoot(), "metrics.txt");
		createBAM(input, SortOrder.unsorted,
				withReadName("indel", Read(0, 1, "10M5I10M"))[0],
				withReadName("sc1", Read(0, 1, "10M5S"))[0],
				withReadName("sc2", Read(0, 1, "10M5S"))[0],
				withReadName("sc2", Read(0, 2, "5S10M5S"))[0],
				withReadName("sc2", Read(0, 2, "5S10M5S"))[0],
				withReadName("indel2", Read(0, 1, "10M5I10M"))[0],
				withReadName("indel2", Read(0, 1, "10M5I10M"))[0]);
		CollectStructuralVariantReadMetrics c = new CollectStructuralVariantReadMetrics();
		c.instanceMain(new String[] {
				"INPUT=" + input.getAbsolutePath(),
				"OUTPUT=" + metricFiles.getAbsolutePath()
		});
		StructuralVariantReadMetrics metrics = getMetrics(metricFiles, StructuralVariantReadMetrics.class);
		assertEquals(2, metrics.INDEL_READS);
		assertEquals(3, metrics.INDEL_READ_ALIGNMENTS);
		assertEquals(2, metrics.SOFT_CLIPPED_READS);
		assertEquals(4, metrics.SOFT_CLIPPED_READ_ALIGNMENTS);
		assertEquals(4, metrics.STRUCTURAL_VARIANT_READS);
		assertEquals(7, metrics.STRUCTURAL_VARIANT_READ_ALIGNMENTS);
	}
}