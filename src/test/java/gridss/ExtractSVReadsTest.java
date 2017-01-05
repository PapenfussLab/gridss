package gridss;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

import au.edu.wehi.idsv.FixedSizeReadPairConcordanceCalculator;
import au.edu.wehi.idsv.IntermediateFilesTest;
import htsjdk.samtools.SAMRecord;

public class ExtractSVReadsTest extends IntermediateFilesTest {
	@Test
	public void hasReadPairingConsistentWithReference_should_search_for_any_concordant_pair() {
		FixedSizeReadPairConcordanceCalculator rpcc = new FixedSizeReadPairConcordanceCalculator(10, 20);
		List<SAMRecord> list = Lists.newArrayList(DP(0, 1, "10M", true, 1, 1, "10M", false));
		assertFalse(ExtractSVReads.hasReadPairingConsistentWithReference(rpcc, list));
		 // positions are consistent but they from alternate mappings of the same read
		list.addAll(Lists.newArrayList(DP(0, 10, "10M", false, 2, 2, "10M", false)));
		assertFalse(ExtractSVReads.hasReadPairingConsistentWithReference(rpcc, list));
		
		list.addAll(Lists.newArrayList(DP(0, 20, "10M", false, 0, 10, "10M", false)));
		assertTrue(ExtractSVReads.hasReadPairingConsistentWithReference(rpcc, list));
	}
	@Test
	public void hasReadAlignmentConsistentWithReference_should_find_any_fully_aligned_mapping_for_segment() {
		List<SAMRecord> list = Lists.newArrayList(DP(0, 1, "10M1S", true, 1, 1, "10M1S", false));
		assertArrayEquals(new boolean[] { false, false }, ExtractSVReads.hasReadAlignmentConsistentWithReference(list));
		list.addAll(Lists.newArrayList(DP(0, 10, "10M", false, 2, 2, "10M1D10M", false)));
		assertArrayEquals(new boolean[] { true, false }, ExtractSVReads.hasReadAlignmentConsistentWithReference(list));
		list.addAll(Lists.newArrayList(DP(0, 10, "5S10M", false, 2, 2, "10X", false)));
		assertArrayEquals(new boolean[] { true, true }, ExtractSVReads.hasReadAlignmentConsistentWithReference(list));
	}
	@Test
	public void hasReadAlignmentConsistentWithReference_should_return_per_segment() {
		List<SAMRecord> list = Lists.newArrayList(withAttr("FI", 2, Read(0, 1, "1M")));
		assertArrayEquals(new boolean[] { false, false, true }, ExtractSVReads.hasReadAlignmentConsistentWithReference(list));
	}
	@Test
	public void should_extract_sv_reads() {
		createInput();
		ExtractSVReads extract = new ExtractSVReads();
		extract.INPUT = input;
		extract.OUTPUT = output;
		extract.METRICS_OUTPUT = new File(output.getAbsolutePath() + ".metrics");
		extract.setup(getHeader(), extract.INPUT);
		extract.acceptFragment(ImmutableList.of(Read(0, 1, "50M50S")), null);
		extract.finish();
		List<SAMRecord> out = getRecords(output);
		assertEquals(1, out.size());
	}
	@Test
	public void should_write_metrics() {
		createInput();
		ExtractSVReads extract = new ExtractSVReads();
		extract.INPUT = input;
		extract.OUTPUT = output;
		extract.METRICS_OUTPUT = new File(output.getAbsolutePath() + ".metrics");
		extract.setup(getHeader(), extract.INPUT);
		extract.acceptFragment(ImmutableList.of(Read(0, 1, "50M50S")), null);
		extract.finish();
		assertTrue(extract.METRICS_OUTPUT.exists());
	}
}
