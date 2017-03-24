package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.junit.Ignore;
import org.junit.Test;
import org.junit.experimental.categories.Category;

import com.google.common.collect.Lists;
import com.google.common.io.Files;

import au.edu.wehi.idsv.alignment.StubFastqAligner;
import au.edu.wehi.idsv.bed.IntervalBed;
import au.edu.wehi.idsv.picard.SynchronousReferenceLookupAdapter;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class SAMEvidenceSourceTest extends IntermediateFilesTest {
	@Test
	public void ensure_metrics_should_write_metrics_files() {
		ProcessingContext pc = getCommandlineContext();
		createInput(RP(0, 100, 200, 100), RP(0, 400, 600, 100));
		SAMEvidenceSource source = new SAMEvidenceSource(pc, input, null, 0);
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
	public void ensure_metrics_should_coverage_blacklist_bed() {
		ProcessingContext pc = getCommandlineContext();
		createInput(RP(0, 100, 200, 100), RP(0, 400, 600, 100));
		SAMEvidenceSource source = new SAMEvidenceSource(pc, input, null, 0);
		source.ensureMetrics();
		assertTrue(pc.getFileSystemContext().getCoverageBlacklistBed(source.getFile()).exists());
	}
	@Test
	public void should_stop_metric_calculation_after_max_records() {
		ProcessingContext pc = getCommandlineContext();
		pc.setCalculateMetricsRecordCount(2);
		createInput(RP(0, 100, 200, 100), RP(0, 400, 600, 100));
		SAMEvidenceSource source = new SAMEvidenceSource(pc, input, null, 0);
		assertEquals(200-100+100, source.getMaxConcordantFragmentSize());
		
		pc.getFileSystemContext().getIdsvMetrics(input).delete();
		pc.getFileSystemContext().getInsertSizeMetrics(input).delete();
		pc.setCalculateMetricsRecordCount(1000);
		source = new SAMEvidenceSource(pc, input, null, 0);
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
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, null, 0);
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
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, null, 0);
		List<DirectedEvidence> list = Lists.newArrayList(source.iterator(new QueryInterval(1, 20, 30)));
		assertTrue(list.stream().allMatch(e ->
			e.getBreakendSummary().overlaps(new BreakendSummary(1, FWD, 20, 20, 30)) ||
			e.getBreakendSummary().overlaps(new BreakendSummary(1, BWD, 20, 20, 30))));
	}
	@Test
	public void should_set_evidence_source_to_self() {
		createInput(Read(0, 1, "50M50S"));
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, null, 0);
		List<DirectedEvidence> list = Lists.newArrayList(source.iterator());
		assertEquals(1, list.size());
		assertEquals(source, list.get(0).getEvidenceSource());
	}
	@Test
	public void should_default_fragment_size_to_read_length_for_unpaired_reads() {
		createInput(Read(0, 1, "50M50S"));
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, null, 0);
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
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, null, 0);
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
				withReadName("r1", Read(1, 1, "15M15S")),
				withReadName("r2", Read(2, 2, "15M15S")),
				withReadName("r3", Read(3, 3, "15M15S")));
		SAMEvidenceSource source = new SAMEvidenceSource(pc, input, null, 0);
		new SplitReadRealigner(pc, new StubFastqAligner(pc)
				.align(inputRecords.get(0), 3, 10, false, "15M")
				.align(inputRecords.get(1), 2, 10, false, "15M")
				.align(inputRecords.get(2), 1, 10, false, "15M"))
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
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, null, 0);
		
		List<DirectedEvidence> result = Lists.newArrayList(source.iterator());
		assertEquals(2, result.size());
	}
	@Test
	public void iterator_should_realign_both_forward_and_backward_soft_clips() throws IOException {
		createInput(
				withReadName("r1", Read(1, 1, "15S15M15S")),
				withReadName("r2", Read(2, 2, "15S15M15S")),
				withReadName("r3", Read(3, 3, "15S15M15S")));
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, null, 0);
		new SplitReadRealigner(getCommandlineContext(), new StubFastqAligner(getCommandlineContext())
				.align(inputRecords.get(0), FWD, 3, 10, false, "15M")
				.align(inputRecords.get(0), BWD, 2, 15, false, "15M")
				.align(inputRecords.get(1), FWD, 2, 10, false, "15M")
				.align(inputRecords.get(2), FWD, 1, 10, false, "15M")
				.align(inputRecords.get(2), BWD, 1, 100, false, "15M"))
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
		SAMEvidenceSource source = new SAMEvidenceSource(pc, input, null, 0);
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
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, null, 0, 0.59);
		
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
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, null, 0, 13, 15);
		
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
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, null, 0, 0, 43);
		List<DirectedEvidence> results = Lists.newArrayList(source.iterator());
		assertEquals(new BreakendSummary(0, BWD, 37, 17, 58), results.get(0).getBreakendSummary());
		assertEquals(new BreakendSummary(0, FWD, 30, 18, 42), results.get(1).getBreakendSummary());
	}
	@Test
	public void iterator_should_sort_sc_with_indel() throws IOException {
		List<SAMRecord> in = new ArrayList<SAMRecord>();
		in.add(Read(0, 2, "5S10M"));
		for (int i = 1; i <= 100; i++) {
			in.add(Read(1, 1, String.format("5M%dD5M5S", i)));
		}
		createInput(in.toArray(new SAMRecord[0]));
		ProcessingContext pc = getCommandlineContext();
		pc.getSoftClipParameters().minLength = 1;
		SAMEvidenceSource source = new SAMEvidenceSource(pc, input, null, 0, 0, 15);
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
		SAMEvidenceSource source = new SAMEvidenceSource(pc, input, null, 0);
		List<DirectedEvidence> result = Lists.newArrayList(source.iterator());
		assertEquals(in.size(), result.size());
		// should be sorted correctly
	}
	@Test
	public void iterator_should_unmap_low_mapq() throws IOException {
		ProcessingContext pc = getCommandlineContext();
		pc.getConfig().minMapq = 10;
		createInput(
				SR(withMapq(10, Read(0, 1, "50M50S"))[0], withMapq(5, Read(1, 1, "50M50S"))[0]).getSAMRecord());
		SAMEvidenceSource source = new SAMEvidenceSource(pc, input, null, 0);
		List<DirectedEvidence> list = Lists.newArrayList(source.iterator());
		assertEquals(1, list.size());
		assertTrue(list.get(0) instanceof SoftClipEvidence);
	}
	@Test
	public void shouldFilter_low_clip_base_qual() {
		SAMRecord r = Read(0, 1, "2M3S");
		r.setBaseQualities(new byte[] { 40, 30, 20, 10, 0 });
		MockSAMEvidenceSource ses = SES();
		ses.getContext().getConfig().getSoftClip().minLength = 1;
		SoftClipEvidence e = SoftClipEvidence.create(ses, BreakendDirection.Forward, r);
		ses.getContext().getConfig().getSoftClip().minAverageQual = 9;
		assertFalse(ses.shouldFilter(e));
		ses.getContext().getConfig().getSoftClip().minAverageQual = 10;
		assertFalse(ses.shouldFilter(e));
		ses.getContext().getConfig().getSoftClip().minAverageQual = 11;
		assertTrue(ses.shouldFilter(e));
	}
	@Test
	@Ignore("base quals now required")
	public void shouldFilter_low_clip_base_qual_should_consider_as_0_if_unknown() {
		SAMRecord r = Read(0, 1, "2M3S");
		r.setBaseQualities(SAMRecord.NULL_QUALS);
		r.setBaseQualities(null);
		MockSAMEvidenceSource ses = SES();
		SoftClipEvidence e = SoftClipEvidence.create(ses, BreakendDirection.Forward, r);
		ses.getContext().getConfig().getSoftClip().minLength = 1;
		ses.getContext().getConfig().getSoftClip().minAverageQual = 0.01f;
		assertTrue(ses.shouldFilter(e));
		r = Read(0, 1, "2M3S");
	}
	@Test
	public void should_filter_dovetailing_reads() {
		SAMEvidenceSource ses = permissiveSES();
		SAMRecord[] rp = RP(0, 100, 100, 20);
		rp[0].setCigarString("10M10S");
		rp[1].setCigarString("10S10M");
		clean(rp[0], rp[1]);
		assertTrue(ses.shouldFilter(rp[0]));
		assertTrue(ses.shouldFilter(rp[1]));
		
		// looks like a dovetail, but is not
		rp = RP(0, 100, 100, 20);
		rp[0].setCigarString("10S10M"); // <-- definitely keep this one
		rp[1].setCigarString("10S10M"); // ideally keep this one too, but we need to know about the mate cigar for that
		clean(rp[0], rp[1]);
		assertFalse(ses.shouldFilter(rp[0]));
		assertFalse(ses.shouldFilter(rp[1]));
	}
	@Test
	public void should_filter_low_complexity_anchors() {
		SAMEvidenceSource ses = permissiveSES();
		ses.getContext().getConfig().minAnchorShannonEntropy = 0.5;
		
		assertTrue(ses.shouldFilter(SoftClipEvidence.create(ses, FWD, withSequence("AT", Read(0, 1, "1M1S"))[0])));
		assertTrue(ses.shouldFilter(SoftClipEvidence.create(ses, FWD, withSequence("AAT", Read(0, 1, "2M1S"))[0])));
		assertTrue(ses.shouldFilter(SoftClipEvidence.create(ses, FWD, withSequence("AAAT", Read(0, 1, "3M1S"))[0])));
		assertTrue(ses.shouldFilter(SoftClipEvidence.create(ses, FWD, withSequence("AAAAT", Read(0, 1, "4M1S"))[0])));
		assertTrue(ses.shouldFilter(SoftClipEvidence.create(ses, FWD, withSequence("AAAAAT", Read(0, 1, "5M1S"))[0])));
		assertTrue(ses.shouldFilter(SoftClipEvidence.create(ses, FWD, withSequence("AAAAAAT", Read(0, 1, "6M1S"))[0])));
		assertTrue(ses.shouldFilter(SoftClipEvidence.create(ses, FWD, withSequence("AAAAAAAT", Read(0, 1, "7M1S"))[0])));
		assertTrue(ses.shouldFilter(SoftClipEvidence.create(ses, FWD, withSequence("AAAAAAAAT", Read(0, 1, "8M1S"))[0])));
		assertTrue(ses.shouldFilter(SoftClipEvidence.create(ses, FWD, withSequence("AAAAAAAAAT", Read(0, 1, "9M1S"))[0])));
		assertTrue(ses.shouldFilter(SoftClipEvidence.create(ses, FWD, withSequence("AAAAAAAAAAT", Read(0, 1, "10M1S"))[0])));
		
		assertTrue(ses.shouldFilter(SoftClipEvidence.create(ses, FWD, withSequence("GT", Read(0, 1, "1M1S"))[0])));
		assertFalse(ses.shouldFilter(SoftClipEvidence.create(ses, FWD, withSequence("GAT", Read(0, 1, "2M1S"))[0])));
		assertFalse(ses.shouldFilter(SoftClipEvidence.create(ses, FWD, withSequence("GAAT", Read(0, 1, "3M1S"))[0])));
		assertFalse(ses.shouldFilter(SoftClipEvidence.create(ses, FWD, withSequence("GAAAT", Read(0, 1, "4M1S"))[0])));
		assertFalse(ses.shouldFilter(SoftClipEvidence.create(ses, FWD, withSequence("GAAAAT", Read(0, 1, "5M1S"))[0])));
		assertFalse(ses.shouldFilter(SoftClipEvidence.create(ses, FWD, withSequence("GAAAAAT", Read(0, 1, "6M1S"))[0])));
		assertFalse(ses.shouldFilter(SoftClipEvidence.create(ses, FWD, withSequence("GAAAAAAT", Read(0, 1, "7M1S"))[0])));
		assertFalse(ses.shouldFilter(SoftClipEvidence.create(ses, FWD, withSequence("GAAAAAAAT", Read(0, 1, "8M1S"))[0])));
		assertFalse(ses.shouldFilter(SoftClipEvidence.create(ses, FWD, withSequence("GAAAAAAAAT", Read(0, 1, "9M1S"))[0])));
		assertTrue(ses.shouldFilter(SoftClipEvidence.create(ses, FWD, withSequence("GAAAAAAAAAT", Read(0, 1, "10M1S"))[0])));
		assertTrue(ses.shouldFilter(SoftClipEvidence.create(ses, FWD, withSequence("GAAAAAAAAAAT", Read(0, 1, "11M1S"))[0])));
		
		assertFalse(ses.shouldFilter(SoftClipEvidence.create(ses, FWD, withSequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTCT", Read(0, 1, "40M1S"))[0]))); // 37	1	1	1	0.503183732
		assertTrue(ses.shouldFilter(SoftClipEvidence.create(ses, FWD, withSequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTCT", Read(0, 1, "41M1S"))[0]))); // 38	1	1	1	0.493619187
		
		assertFalse(ses.shouldFilter(SoftClipEvidence.create(ses, FWD, withSequence("TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTCT", Read(0, 1, "1S40M1S"))[0]))); // 37	1	1	1	0.503183732
		assertTrue(ses.shouldFilter(SoftClipEvidence.create(ses, FWD, withSequence("TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTCT", Read(0, 1, "1S41M1S"))[0]))); // 38	1	1	1	0.493619187
	}
	@Test
	public void should_filter_reference_supporting_reads() {
		assertTrue(SES().shouldFilter(SR(Read(0, 10, "5M5S"), Read(1, 10, "5M"))));
		assertTrue(SES().shouldFilter(SR(Read(1, 10, "5M5S"), Read(0, 10, "5M"))));
	}
	@Test
	public void should_ignore_duplicates() {
		SAMRecord r = Read(1, 1, "10M5S");
		assertFalse(SES().shouldFilter(r)); // precondition
		r.setDuplicateReadFlag(true);
		assertTrue(SES().shouldFilter(r));
	}
	@Test
	public void blacklisted_alignment_should_be_unmapped() {
		ProcessingContext pc = getContext();
		MockSAMEvidenceSource ses = SES(pc);
		IntervalBed blacklist = new IntervalBed(pc.getDictionary(), pc.getLinear());
		blacklist.addInterval(0, 10, 12);
		ses.setBlacklistedRegions(blacklist);
		assertFalse(ses.transform(Read(1, 10, "1M")).getReadUnmappedFlag());
		assertFalse(ses.transform(Read(0, 9, "1M")).getReadUnmappedFlag());
		assertTrue(ses.transform(Read(0, 10, "1M")).getReadUnmappedFlag());
		assertTrue(ses.transform(Read(0, 11, "1M")).getReadUnmappedFlag());
		assertTrue(ses.transform(Read(0, 12, "1M")).getReadUnmappedFlag());
		assertFalse(ses.transform(Read(0, 13, "1M")).getReadUnmappedFlag());
	}
	@Test
	public void blacklist_mate_should_be_treated_as_unmapped() {
		ProcessingContext pc = getContext();
		MockSAMEvidenceSource ses = SES(pc);
		IntervalBed blacklist = new IntervalBed(pc.getDictionary(), pc.getLinear());
		blacklist.addInterval(0, 10, 12);
		ses.setBlacklistedRegions(blacklist);
		assertFalse(ses.transform(DP(0, 1, "1M", true, 1, 10, "1M", false)[0]).getMateUnmappedFlag());
		assertFalse(ses.transform(DP(0, 1, "1M", true, 0, 9, "1M", false)[0]).getMateUnmappedFlag());
		assertTrue(ses.transform(DP(0, 1, "1M", true, 0, 10, "1M", false)[0]).getMateUnmappedFlag());
		assertTrue(ses.transform(DP(0, 1, "1M", true, 0, 11, "1M", false)[0]).getMateUnmappedFlag());
		assertTrue(ses.transform(DP(0, 1, "1M", true, 0, 12, "1M", false)[0]).getMateUnmappedFlag());
		assertFalse(ses.transform(DP(0, 1, "1M", true, 0, 13, "1M", false)[0]).getMateUnmappedFlag());
	}
	@Test
	public void blacklist_chimeric_alignment_should_be_treated_as_unmapped() {
		ProcessingContext pc = getContext();
		MockSAMEvidenceSource ses = SES(pc);
		IntervalBed blacklist = new IntervalBed(pc.getDictionary(), pc.getLinear());
		blacklist.addInterval(0, 10, 12);
		ses.setBlacklistedRegions(blacklist);
		SAMRecord read = Read(0, 1, "1M1S");
		read.setAttribute("SA", "polyA,10,+,1S1M,10,10");
		assertEquals("", ses.transform(read).getAttribute("SA"));
	}
	@Test
	public void should_filter_breakpoints_touching_blacklist() {
		ProcessingContext pc = getContext();
		MockSAMEvidenceSource ses = SES(pc);
		IntervalBed blacklist = new IntervalBed(pc.getDictionary(), pc.getLinear());
		blacklist.addInterval(0, 10, 12);
		ses.setBlacklistedRegions(blacklist);
		assertTrue(ses.shouldFilter(NRRP(ses, DP(0, 1, "1M", true, 1, 100, "1M", false))));
		assertTrue(ses.shouldFilter(NRRP(ses, DP(1, 100, "1M", false, 0, 1, "1M", true))));
	}
	//@Test
	@Category(Hg38Tests.class)
	public void regression_indel_should_not_throw_exception() throws IOException {
		File ref = Hg38Tests.findHg38Reference();
		ProcessingContext pc = new ProcessingContext(getFSContext(), ref, new SynchronousReferenceLookupAdapter(new IndexedFastaSequenceFile(ref)), null, getConfig());
		Files.copy(new File("src/test/resources/indelerror.bam"), input);
		SAMEvidenceSource source = new SAMEvidenceSource(pc, input, null, 0);
		List<DirectedEvidence> reads = Lists.newArrayList(source.iterator());
	}
	@Test
	public void should_filter_reads_aligning_outside_contig_bounds() {
		SAMEvidenceSource ses = permissiveSES();
		assertFalse(ses.shouldFilter(Read(0, 1, "10M")));
		assertFalse(ses.shouldFilter(Read(0, 9999, "2M")));
		assertTrue(ses.shouldFilter(Read(0, 0, "10M")));
		assertTrue(ses.shouldFilter(Read(0, 10000, "2M")));
	}
	@Test
	public void should_not_filter_XNX_placholder_outside_of_bounds() {
		SAMEvidenceSource ses = permissiveSES();
		assertFalse(ses.shouldFilter(Read(0, -10, "10S1X1N1X")));
	}
}
