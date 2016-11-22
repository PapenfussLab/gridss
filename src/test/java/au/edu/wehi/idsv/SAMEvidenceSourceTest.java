package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.Lists;

import au.edu.wehi.idsv.alignment.StubFastqAligner;
import au.edu.wehi.idsv.bed.IntervalBed;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;

public class SAMEvidenceSourceTest extends IntermediateFilesTest {
	@Test
	public void ensure_metrics_should_write_metrics_files() {
		ProcessingContext pc = getCommandlineContext();
		createInput(RP(0, 100, 200, 100), RP(0, 400, 600, 100));
		SAMEvidenceSource source = new SAMEvidenceSource(pc, input, 0);
		source.ensureMetrics();
		assertTrue(pc.getFileSystemContext().getIdsvMetrics(source.getFile()).exists());
		assertTrue(pc.getFileSystemContext().getInsertSizeMetrics(source.getFile()).exists());
		assertTrue(pc.getFileSystemContext().getCigarMetrics(source.getFile()).exists());
		assertTrue(pc.getFileSystemContext().getMapqMetrics(source.getFile()).exists());
		assertNotNull(source.getMetrics());
		assertNotNull(source.getMetrics().getCigarDetailMetrics());
		assertNotNull(source.getMetrics().getCigarDistribution());
		assertNotNull(source.getMetrics().getIdsvMetrics());
		assertNotNull(source.getMetrics().getInsertSizeDistribution());
		assertNotNull(source.getMetrics().getMapqMetrics());
	}
	@Test
	public void should_stop_metric_calculation_after_max_records() {
		ProcessingContext pc = getCommandlineContext();
		pc.setCalculateMetricsRecordCount(2);
		createInput(RP(0, 100, 200, 100), RP(0, 400, 600, 100));
		SAMEvidenceSource source = new SAMEvidenceSource(pc, input, 0);
		assertEquals(200-100+100, source.getMaxConcordantFragmentSize());
		
		pc.getFileSystemContext().getIdsvMetrics(input).delete();
		pc.getFileSystemContext().getInsertSizeMetrics(input).delete();
		pc.setCalculateMetricsRecordCount(1000);
		source = new SAMEvidenceSource(pc, input, 0);
		assertEquals(600-400+100, source.getMaxConcordantFragmentSize());
	}
	@Test
	public void iterator_should_return_all_evidence() {
		createInput(new SAMRecord[] { Read(1, 1, "50M50S") },
				RP(0, 100, 200, 100),
				DP(1, 1, "100M", true, 2, 5, "100M", true),
			   DP(1, 2, "100M", true, 2, 4, "100M", true),
			   DP(1, 3, "100M", true, 2, 6, "100M", true),
			   OEA(1, 4, "100M", false));
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, 0);
		List<DirectedEvidence> list = Lists.newArrayList(source.iterator());
		assertEquals(8, list.size()); // 1 SC + 3 * 2 DP + 1 OEA
	}
	@Test
	public void iterator_filter_to_breakends_overlapping_query_interval() {
		List<SAMRecord> in = new ArrayList<>();
		for (int i = 1; i < 100; i++) {
			in.add(Read(1, i, "5S5M"));
			in.add(Read(1, i, "5M2I5M"));
			in.add(Read(1, i, "5M5S"));
			in.add(Read(1, i, "1X2N1X5S"));
			in.add(Read(1, i, "5S1X2N1X"));
			Collections.addAll(in, RP(0, i, i + 10, 5));
			Collections.addAll(in, RP(1, i, i + 10, 5));
			Collections.addAll(in, OEA(1, i, "5M", true));
			Collections.addAll(in, OEA(1, i, "5M", false));
			Collections.addAll(in, DP(1, i, "5M", true, 0, 1, "5M", false));
			Collections.addAll(in, DP(1, i, "5M", false, 0, 1, "5M", false));
		}
		createInput(in);
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, 0);
		List<DirectedEvidence> list = Lists.newArrayList(source.iterator(new QueryInterval(1, 20, 30)));
		assertTrue(list.stream().allMatch(e ->
			e.getBreakendSummary().overlaps(new BreakendSummary(1, FWD, 20, 20, 30)) ||
			e.getBreakendSummary().overlaps(new BreakendSummary(1, BWD, 20, 20, 30))));
	}
	@Test
	public void should_set_evidence_source_to_self() {
		createInput(Read(0, 1, "50M50S"));
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, 0);
		List<DirectedEvidence> list = Lists.newArrayList(source.iterator());
		assertEquals(1, list.size());
		assertEquals(source, list.get(0).getEvidenceSource());
	}
	@Test
	public void should_default_fragment_size_to_read_length_for_unpaired_reads() {
		createInput(Read(0, 1, "50M50S"));
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, 0);
		assertEquals(100, source.getMaxConcordantFragmentSize());
	}
	@Test
	public void iterator_evidence_should_be_sorted_by_evidence_natural_ordering() {
		createInput(new SAMRecord[] { Read(1, 1, "50M50S") },
				new SAMRecord[] { Read(1, 1, "50S50M") },
				RP(0, 100, 200, 100),
				DP(1, 1, "100M", true, 2, 5, "100M", true),
			   DP(1, 2, "100M", false, 2, 4, "100M", true),
			   DP(1, 3, "100M", true, 2, 6, "100M", true),
			   OEA(1, 4, "100M", false),
			   OEA(1, 4, "100M", true));
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, 0);
		List<DirectedEvidence> result = Lists.newArrayList(source.iterator());
		List<DirectedEvidence> sorted = Lists.newArrayList(result);
		Collections.sort(sorted, DirectedEvidenceOrder.ByNatural);
		assertEquals(sorted, result);
	}
	@Test
	public void iterator_should_return_both_sides_of_split_reads() throws IOException {
		ProcessingContext pc = getCommandlineContext();
		pc.getConfig().minMapq = 10;
		createInput(
				withReadName("r1", Read(0, 1, "15M15S")),
				withReadName("r2", Read(1, 2, "15M15S")),
				withReadName("r3", Read(2, 3, "15M15S")));
		SAMEvidenceSource source = new SAMEvidenceSource(pc, input, 0);
		new SplitReadRealigner(pc, new StubFastqAligner(pc)
				.align(inputRecords.get(0), 2, 10, false, "15M")
				.align(inputRecords.get(1), 1, 10, false, "15M")
				.align(inputRecords.get(2), 0, 10, false, "15M"))
			.createSupplementaryAlignments(input, input);
		
		List<DirectedEvidence> result = Lists.newArrayList(source.iterator());
		assertEquals(6, result.size());
		assertTrue(result.get(0) instanceof SplitReadEvidence);
		assertTrue(result.get(1) instanceof SplitReadEvidence);
		assertTrue(result.get(2) instanceof SplitReadEvidence);
	}
	@Test
	public void iterator_should_iterator_over_both_forward_and_backward_soft_clips() {
		createInput(withReadName("r1", Read(0, 1, "15S15M15S")));
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, 0);
		
		List<DirectedEvidence> result = Lists.newArrayList(source.iterator());
		assertEquals(2, result.size());
	}
	@Test
	public void iterator_should_realign_both_forward_and_backward_soft_clips() throws IOException {
		createInput(
				withReadName("r1", Read(0, 1, "15S15M15S")),
				withReadName("r2", Read(1, 2, "15S15M15S")),
				withReadName("r3", Read(2, 3, "15S15M15S")));
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, 0);
		new SplitReadRealigner(getCommandlineContext(), new StubFastqAligner(getCommandlineContext())
				.align(inputRecords.get(0), FWD, 2, 10, false, "15M")
				.align(inputRecords.get(0), BWD, 1, 15, false, "15M")
				.align(inputRecords.get(1), FWD, 1, 10, false, "15M")
				.align(inputRecords.get(2), FWD, 0, 10, false, "15M")
				.align(inputRecords.get(2), BWD, 0, 100, false, "15M"))
			.createSupplementaryAlignments(input, input);
		
		List<DirectedEvidence> result = Lists.newArrayList(source.iterator());
		assertEquals(6 + 5, result.size());
		// 4 split reads = 8 breakends
		assertEquals(5 * 2, result.stream().filter(e -> e instanceof SplitReadEvidence).count());
		// and a soft clipped read
		assertEquals(1, result.stream().filter(e -> e instanceof SoftClipEvidence).count());
	}
	@Test
	public void remote_soft_clips_should_match_location() throws IOException {
		ProcessingContext pc = getCommandlineContext();
		pc.getConfig().minMapq = 10;
		createInput(
				withReadName("r2", Read(1, 2, "15M15S")));
		SAMEvidenceSource source = new SAMEvidenceSource(pc, input, 0);
		new SplitReadRealigner(getCommandlineContext(), new StubFastqAligner(getCommandlineContext())
				.align(inputRecords.get(0), FWD, 1, 10, false, "15M"))
			.createSupplementaryAlignments(input, input); 
		
		List<DirectedEvidence> result = Lists.newArrayList(source.iterator());
		assertEquals(result.get(0).getBreakendSummary(), ((BreakpointSummary)result.get(1).getBreakendSummary()).remoteBreakpoint());
	}
	@Test
	public void should_construct_percentage_based_calculator() {
		createInput(
				RP(0, 100, 110, 5), // 15
				RP(0, 100, 111, 5), 
				RP(0, 100, 112, 5), // --- 17
				RP(0, 100, 113, 5), // ---
				RP(0, 100, 114, 5), // --- concordant 
				RP(0, 100, 115, 5), // --- reads
				RP(0, 100, 116, 5), // ---
				RP(0, 100, 117, 5), // --- 22
				RP(0, 100, 118, 5),
				RP(0, 100, 119, 5) // 24
				);
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, 0, 0.59);
		
		assertEquals(PercentageReadPairConcordanceCalculator.class, source.getReadPairConcordanceCalculator().getClass());
		assertEquals(22, source.getMaxConcordantFragmentSize());
		assertEquals(5, source.getMaxReadLength());
		
		List<DirectedEvidence> list = Lists.newArrayList(source.iterator());
		assertEquals(4 * 2, list.size()); // 4 discordant pairs
	}
	@Test
	public void should_construct_fixed_calculator() {
		createInput(
				RP(0, 1, 12, 1), // 12
				RP(0, 1, 13, 1), // 13 conc
				RP(0, 1, 14, 1), // 14 conc
				RP(0, 1, 15, 1), // 15 conc
				RP(0, 1, 16, 1)); // 16
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, 0, 13, 15);
		
		assertEquals(FixedSizeReadPairConcordanceCalculator.class, source.getReadPairConcordanceCalculator().getClass());
		assertEquals(15, source.getMaxConcordantFragmentSize());
		assertEquals(1, source.getMaxReadLength());
		List<DirectedEvidence> list = Lists.newArrayList(source.iterator());
		assertEquals(2 * 2, list.size());
	}

	@Test
	public void iterator_should_sort_worst_case_read_pairs() {
		//          1         2         3         4         5         6         7
		// 1234567890123456789012345678901234567890123456789012345678901234567890
		// MMMMMMMMMMMMMMMMMM>
		// ^--read length --^
		// ^--------max concordant fragment size-----^
		//                  |-----breakend call-----|
		//                 |---------------breakend call------------|
		//                                                <SSSSSSSSSM
		//                ^--------max concordant fragment size-----^
		// ^ alignment start                                          ^ alignment start
		List<SAMRecord> in = new ArrayList<SAMRecord>();
		for (int i = 1; i < 58; i++) {
			in.addAll(Lists.newArrayList(OEA(0, i, "18M", true)));
		}
		in.addAll(Lists.newArrayList(OEA(0, 58, "17S1M", false)));
		createInput(in.toArray(new SAMRecord[0]));
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, 0, 0, 43);
		List<DirectedEvidence> results = Lists.newArrayList(source.iterator());
		assertEquals(new BreakendSummary(0, BWD, 37, 17, 58), results.get(0).getBreakendSummary());
		assertEquals(new BreakendSummary(0, FWD, 30, 18, 42), results.get(1).getBreakendSummary());
	}
	@Test
	public void iterator_should_sort_sc_with_indel() throws IOException {
		List<SAMRecord> in = new ArrayList<SAMRecord>();
		in.add(Read(0, 2, "5S10M"));
		for (int i = 1; i <= 100; i++) {
			in.add(Read(0, 1, String.format("5M%dD5M5S", i)));
		}
		createInput(in.toArray(new SAMRecord[0]));
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, 0, 0, 15);
		List<DirectedEvidence> results = Lists.newArrayList(source.iterator());
		for (int i = 1; i < results.size(); i++) {
			assertTrue(results.get(i-1).getBreakendSummary().start <= results.get(i).getBreakendSummary().start);
		}
		assertEquals(3 * 100 + 1, results.size());
	}
	@Test
	public void sort_window_should_allow_for_microhomology_causing_bounds_change() {
		ProcessingContext pc = getCommandlineContext();
		pc.getConfig().minMapq = 0;
		pc.getSoftClipParameters().minLength = 1;
		pc.getSoftClipParameters().minAnchorIdentity = 0;
		pc.getConfig().adapters = new AdapterHelper(new String[0]);
		List<SAMRecord> in = Lists.newArrayList();
		List<SAMRecord> realn = Lists.newArrayList();
		for (int i = 100; i < 200; i++) {
			in.add(withReadName("a_" + Integer.toString(i), withSequence("NNNNNNNNNNN", Read(0, i, "10M1S")))[0]);
			in.add(withReadName("b_" + Integer.toString(i), Read(0, i, "10S1M"))[0]);
			in.add(withReadName("c_" + Integer.toString(i), withSequence("NNNNNNNNNNN", Read(0, i, "1S10M")))[0]);
			in.add(withReadName("d_" + Integer.toString(i), withSequence("NNNNNNNNNNN", Read(0, i, "10M1S")))[0]);
			realn.add(withReadName(String.format("0#%d#0#fa_%d", i, i), withSequence("N", Read(1, 10, "1M")))[0]);
			realn.add(withReadName(String.format("0#%d#0#bb_%d", i, i), withSequence("NNNNNNNNNN", Read(1, 10, "10M")))[0]);
			realn.add(withReadName(String.format("0#%d#0#bc_%d", i, i), withSequence("N", Read(1, 10, "1M")))[0]);
			realn.add(withReadName(String.format("0#%d#0#fd_%d", i, i), withSequence("NNNNNNNNNN", Read(1, 10, "10M")))[0]);
		}
		createInput(in);
		SAMEvidenceSource source = new SAMEvidenceSource(pc, input, 0);
		List<DirectedEvidence> result = Lists.newArrayList(source.iterator());
		assertEquals(in.size(), result.size());
		// should be sorted correctly
	}
	@Test
	public void should_filter_blacklisted_regions() {
		ProcessingContext pc = getCommandlineContext();
		IntervalBed blacklist = new IntervalBed(pc.getDictionary(), pc.getLinear());
		blacklist.addInterval(2, 1, 1000);
		pc.setBlacklistedRegions(blacklist);
		createInput(
				new SAMRecord[] {Read(1, 1, "50M50S") },
				RP(0, 200, 100), // max frag size
				DP(1, 1, "100M", true, 2, 5, "100M", true),
			   DP(1, 2, "100M", true, 2, 4, "100M", true),
			   DP(1, 3, "100M", true, 2, 6, "100M", true),
			   OEA(1, 4, "100M", false));
		SAMEvidenceSource source = new SAMEvidenceSource(pc, input, 0);
		List<DirectedEvidence> list = Lists.newArrayList(source.iterator());
		assertEquals(2, list.size()); // SC & OEA not on (2)
	}
}
