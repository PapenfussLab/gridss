package gridss;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import org.junit.Assert;
import org.junit.Test;
import org.junit.experimental.categories.Category;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.io.Files;

import au.edu.wehi.idsv.FixedSizeReadPairConcordanceCalculator;
import au.edu.wehi.idsv.Hg38Tests;
import au.edu.wehi.idsv.IntermediateFilesTest;
import au.edu.wehi.idsv.ReadPairConcordanceMethod;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.picard.SynchronousReferenceLookupAdapter;
import au.edu.wehi.idsv.sam.ChimericAlignment;
import gridss.analysis.StructuralVariantReadMetrics;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

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
		StructuralVariantReadMetrics metric = Iterators.getOnlyElement(Iterables.filter(MetricsFile.readBeans(extract.METRICS_OUTPUT), StructuralVariantReadMetrics.class).iterator(), null);
		assertEquals(1, metric.SOFT_CLIPPED_READS);
	}
	/*
	@Test
	public void should_not_extract_unclipped_alignment_overlapping_blacklist() {
		createInput();
		IntervalBed blacklist = new IntervalBed(getSequenceDictionary(), new PaddedLinearGenomicCoordinate(getSequenceDictionary()));
		blacklist.addInterval(0, 2, 2);
				
		ExtractSVReads extract = new ExtractSVReads();
		extract.INPUT = input;
		extract.OUTPUT = output;
		extract.setBlacklist(blacklist);
		extract.setup(getHeader(), extract.INPUT);
		extract.acceptFragment(ImmutableList.of(Read(0, 1, "1M1S")), null); // Soft clip breakpoint position overlaps
		extract.acceptFragment(ImmutableList.of(Read(0, 2, "1M1S")), null); // overlaps
		extract.acceptFragment(ImmutableList.of(Read(0, 3, "1M1S")), null);
		extract.acceptFragment(ImmutableList.of(Read(0, 4, "1M1S")), null);
		extract.acceptFragment(ImmutableList.of(Read(0, 3, "1S1M")), null); // overlaps
		extract.finish();
		List<SAMRecord> out = getRecords(output);
		assertEquals(2, out.size());
	}
	@Test
	public void should_not_extract_any_read_from_fragment_overlapping_read_pair_breakend() {
		createInput();
		IntervalBed blacklist = new IntervalBed(getSequenceDictionary(), new PaddedLinearGenomicCoordinate(getSequenceDictionary()));
		blacklist.addInterval(0, 2, 2);
				
		ExtractSVReads extract = new ExtractSVReads();
		extract.INPUT = input;
		extract.OUTPUT = output;
		extract.FIXED_READ_PAIR_CONCORDANCE_MAX_FRAGMENT_SIZE = 2;
		extract.FIXED_READ_PAIR_CONCORDANCE_MIN_FRAGMENT_SIZE = 2;
		extract.READ_PAIR_CONCORDANCE_METHOD = ReadPairConcordanceMethod.FIXED;
		extract.setBlacklist(blacklist);
		extract.setup(getHeader(), extract.INPUT);
		extract.acceptFragment(Lists.newArrayList(DP(0, 1, "1M", true, 1, 1, "1M", false)), null); // breakend overlaps
		extract.acceptFragment(Lists.newArrayList(DP(0, 1, "1M", false, 1, 1, "1M", false)), null); // other direction = no overlap 
		extract.acceptFragment(Lists.newArrayList(DP(0, 2, "1M", true, 1, 1, "1M", false)), null); // read overlaps
		extract.acceptFragment(Lists.newArrayList(DP(0, 2, "1M", false, 1, 1, "1M", false)), null); // read overlaps 
		extract.acceptFragment(Lists.newArrayList(DP(0, 3, "1M", true, 1, 1, "1M", false)), null); // other direction = no overlap 
		extract.acceptFragment(Lists.newArrayList(DP(0, 3, "1M", false, 1, 1, "1M", false)), null); // breakend overlaps
		extract.finish();
		List<SAMRecord> out = getRecords(output);
		assertEquals(2 * 2, out.size());
	}
	@Test
	public void should_not_extract_any_split_reads_overlapping_blacklist() {
		fail();
	}
	*/
	@Test
	public void hasReadAlignmentConsistentWithReference_should_consider_SA_alignments() {
		SAMRecord r = Read(0, 1, "5M5S");
		r.setAttribute("SA", new ChimericAlignment(Read(0, 1, "10M")).toString());
		Assert.assertTrue(ExtractSVReads.hasReadAlignmentConsistentWithReference(ImmutableList.of(r))[0]);
	}
	@Test
	public void hasReadPairingConsistentWithReference_should_use_primary_alignment_for_split_alignments() {
		SAMRecord r = Read(0, 1, "5M5S");
		r.setAttribute("SA", new ChimericAlignment(Read(1, 50, "10M")).toString());
		r.setMateAlignmentStart(100);
		r.setMateNegativeStrandFlag(true);
		r.setMateReferenceIndex(1);
		r.setReadPairedFlag(true);
		r.setMateUnmappedFlag(false);
		r.setSupplementaryAlignmentFlag(true);
		Assert.assertTrue(ExtractSVReads.hasReadPairingConsistentWithReference(new FixedSizeReadPairConcordanceCalculator(0, 100), ImmutableList.of(r)));
	}
	@Test
	public void should_extract_fully_mapped_split_read() {
		ExtractSVReads extract = new ExtractSVReads();
		extract.instanceMain(new String[] {
				"INPUT=src/test/resources/fullymappedsplitread.sam",
				"OUTPUT=" + output.getAbsolutePath()
		});
		List<SAMRecord> records = getRecords(output);
		assertEquals(0, records.size());
	}
	@Test
	@Category(Hg38Tests.class)
	public void regression_should_extract_split_read_alignments_as_group() throws IOException {
		File ref = Hg38Tests.findHg38Reference();
		ReferenceLookup lookup = new SynchronousReferenceLookupAdapter(new IndexedFastaSequenceFile(ref));
		Files.copy(new File("src/test/resources/sa.split/test1.sam"), input);
		
		ExtractSVReads extract = new ExtractSVReads();
		extract.setReference(lookup); //new ProcessingContext(getFSContext(), ref, lookup, null, getConfig());
		extract.MIN_CLIP_LENGTH = 4;
		extract.INSERT_SIZE_METRICS = new File("src/test/resources/sa.split/test.insert_size_metrics");
		extract.READ_PAIR_CONCORDANCE_METHOD = ReadPairConcordanceMethod.PERCENTAGE;
		extract.OUTPUT = output;
		extract.INPUT = input;
		try (SamReader reader = SamReaderFactory.make().open(input)) {
			extract.setup(reader.getFileHeader(), input);
		}
		List<SAMRecord> records = getRecords(input);
		List<Boolean> result = records.stream().map(r -> extract.shouldExtract(ImmutableList.of(r), lookup)[0]).collect(Collectors.toList());
		/*
		for (int i = 0; i < records.size(); i++) {
			SAMRecord r = records.get(i);
			System.out.print(r.getSupplementaryAlignmentFlag() ? "S" : " ");
			System.out.print(r.getFirstOfPairFlag() ? "1" : "2");
			System.out.print(" Extracted=" + (result.get(i) ? "y" : "n"));
			System.out.print(" HasConcPair=" + (ExtractSVReads.hasReadPairingConsistentWithReference(extract.getReadPairConcordanceCalculator(), ImmutableList.of(r)) ? "y" : "n"));
			boolean[] ra = ExtractSVReads.hasReadAlignmentConsistentWithReference(ImmutableList.of(r));
			System.out.print(" HasConcRead=" + (ra[0] ? "y" : "n") + (ra[1] ? "y" : "n"));
			System.out.println(" " + new ChimericAlignment(r).toString());
		}
		*/
		Assert.assertEquals(records.stream().map(r -> r.getSecondOfPairFlag()).collect(Collectors.toList()), result);
		lookup.close();
	}
	@Test
	@Category(Hg38Tests.class)
	public void regression_should_extract_split_read_alignments_as_group_2() throws IOException {
		File ref = Hg38Tests.findHg38Reference();
		ReferenceLookup lookup = new SynchronousReferenceLookupAdapter(new IndexedFastaSequenceFile(ref));
		Files.copy(new File("src/test/resources/sa.split/test2.sam"), input);
		
		ExtractSVReads extract = new ExtractSVReads();
		extract.setReference(lookup); //new ProcessingContext(getFSContext(), ref, lookup, null, getConfig());
		extract.MIN_CLIP_LENGTH = 4;
		extract.INSERT_SIZE_METRICS = new File("src/test/resources/sa.split/test.insert_size_metrics");
		extract.READ_PAIR_CONCORDANCE_METHOD = ReadPairConcordanceMethod.PERCENTAGE;
		extract.OUTPUT = output;
		extract.INPUT = input;
		try (SamReader reader = SamReaderFactory.make().open(input)) {
			extract.setup(reader.getFileHeader(), input);
		}
		List<SAMRecord> records = getRecords(input);
		List<Boolean> result = records.stream().map(r -> extract.shouldExtract(ImmutableList.of(r), lookup)[0]).collect(Collectors.toList());
		/*
		for (int i = 0; i < records.size(); i++) {
			SAMRecord r = records.get(i);
			System.out.print(r.getSupplementaryAlignmentFlag() ? "S" : " ");
			System.out.print(r.getFirstOfPairFlag() ? "1" : "2");
			System.out.print(" Extracted=" + (result.get(i) ? "y" : "n"));
			System.out.print(" HasConcPair=" + (ExtractSVReads.hasReadPairingConsistentWithReference(extract.getReadPairConcordanceCalculator(), ImmutableList.of(r)) ? "y" : "n"));
			boolean[] ra = ExtractSVReads.hasReadAlignmentConsistentWithReference(ImmutableList.of(r));
			System.out.print(" HasConcRead=" + (ra[0] ? "y" : "n") + (ra[1] ? "y" : "n"));
			System.out.println(" " + new ChimericAlignment(r).toString());
		}*/
		// primary read pair alignment implies an unexpected library fragment size = extract them all  
		Assert.assertEquals(ImmutableList.of(true, true, true, true), result);
		lookup.close();
	}
	@Test
	public void should_extract_sam_flag_without_insert_size_metrics() {
		ExtractSVReads extract = new ExtractSVReads();
		extract.instanceMain(new String[] {
			"INPUT=example/chr12.1527326.DEL1024.bam",
			"OUTPUT=" + output.toString(),
			"READ_PAIR_CONCORDANCE_METHOD=SAM_FLAG",
			"METRICS_OUTPUT=" + output.toString() + ".sv.metrics",
		});
		assertTrue(output.exists());
	}
}
