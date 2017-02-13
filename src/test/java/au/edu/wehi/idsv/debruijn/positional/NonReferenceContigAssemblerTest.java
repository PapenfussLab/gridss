package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

import org.junit.Assert;
import org.junit.Test;

import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;

import au.edu.wehi.idsv.AssemblyAttributes;
import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.DiscordantReadPair;
import au.edu.wehi.idsv.NonReferenceReadPair;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.SequentialIdGenerator;
import au.edu.wehi.idsv.SoftClipEvidence;
import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.SAMRecord;


public class NonReferenceContigAssemblerTest extends TestHelper {
	private NonReferenceContigAssembler caller;
	private EvidenceTracker tracker;
	public List<SAMRecord> go(ProcessingContext pc, boolean collapse, DirectedEvidence... input) {
		caller = create(pc, collapse, input);
		List<SAMRecord> assemblies = Lists.newArrayList(caller);
		return assemblies;
	}
	private NonReferenceContigAssembler create(ProcessingContext pc, boolean collapse, DirectedEvidence... input) {
		Arrays.sort(input, DirectedEvidence.ByStartEnd);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, Arrays.stream(input).map(e -> (SAMEvidenceSource)e.getEvidenceSource()).collect(Collectors.toList()), null);
		int maxEvidenceWidth = aes.getMaxConcordantFragmentSize() - aes.getMinConcordantFragmentSize() + 1;
		int maxReadLength = maxReadLength(input);
		int k = pc.getAssemblyParameters().k;
		int maxPathLength = pc.getAssemblyParameters().positional.maxPathLengthInBases(maxReadLength);
		int maxPathCollapseLength = pc.getAssemblyParameters().errorCorrection.maxPathCollapseLengthInBases(maxReadLength);
		tracker = new EvidenceTracker();
		SupportNodeIterator supportIt = new SupportNodeIterator(k, Arrays.stream(input).iterator(), aes.getMaxConcordantFragmentSize(), tracker, pc.getAssemblyParameters().includePairAnchors, pc.getAssemblyParameters().pairAnchorMismatchIgnoreEndBases);
		AggregateNodeIterator agIt = new AggregateNodeIterator(supportIt);
		Iterator<KmerPathNode> pnIt = new PathNodeIterator(agIt, maxPathLength, k);
		if (collapse) {
			pnIt = new PathCollapseIterator(pnIt, k, maxPathCollapseLength, pc.getAssemblyParameters().errorCorrection.maxBaseMismatchForCollapse, pc.getAssemblyParameters().errorCorrection.collapseBubblesOnly, 0);
			pnIt = new PathSimplificationIterator(pnIt, maxPathLength, maxEvidenceWidth);
		}
		caller = new NonReferenceContigAssembler(pnIt, 0, maxEvidenceWidth + maxReadLength + 2, maxReadLength, k, aes, new SequentialIdGenerator("asm"), tracker, "test");
		return caller;
	}
	@Test
	public void should_call_simple_fwd_SC() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SoftClipEvidence sce = SCE(FWD, withSequence("ACGTGGTCGACC", Read(0, 5, "6M6S")));
		List<SAMRecord> output = go(pc, true, sce);
		assertEquals(1, output.size());
		assertEquals("ACGTGGTCGACC", S(output.get(0).getReadBases()));
		assertEquals("6M6S", output.get(0).getCigarString());
		assertEquals(0, (int)output.get(0).getReferenceIndex());
		assertEquals(5, output.get(0).getAlignmentStart());
	}
	@Test
	public void should_set_assembly_attributes() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SoftClipEvidence sce = SCE(FWD, withSequence("ACGTGGTCGACC", Read(0, 5, "6M6S")));
		List<SAMRecord> output = go(pc, true, sce);
		AssemblyAttributes attr = new AssemblyAttributes(output.get(0));
		assertNotNull(attr);
		assertEquals(1, attr.getAssemblySupportCount());
	}
	@Test
	public void should_call_simple_bwd_SC() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SoftClipEvidence sce = SCE(BWD, withSequence("ACGTGGTCGACC", Read(0, 10, "6S6M")));
		List<SAMRecord> output = go(pc, true, sce);
		
		assertEquals(1, output.size());
		assertEquals("ACGTGGTCGACC", S(output.get(0).getReadBases()));
		assertEquals("6S6M", output.get(0).getCigarString());
		assertEquals(0, (int)output.get(0).getReferenceIndex());
		assertEquals(10, output.get(0).getAlignmentStart());
	}
	@Test
	public void should_call_overlapping_SC() {
		// 12345678901234567890
		// ACGTGGCTCGACC
		// ------
		//        ------
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SoftClipEvidence fwd = SCE(BWD, withSequence("ACGTGGCTCGACC", Read(0, 8, "7S6M")));
		SoftClipEvidence bwd = SCE(FWD, withSequence("ACGTGGCTCGACC", Read(0, 1, "6M7S")));
		List<SAMRecord> output = go(pc, true, fwd, bwd);
		assertEquals(1, output.size());
		assertEquals("6M1I1D6M", output.get(0).getCigarString());
	}
	@Test
	public void should_call_simple_fwd_RP() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SAMEvidenceSource ses = SES(100, 200);
		DiscordantReadPair e = (DiscordantReadPair)NRRP(ses, withSequence("ACGTGGTCGACC", DP(0, 50, "12M", true, 1, 1, "12M", false)));
		List<SAMRecord> output = go(pc, true, e);
		assertEquals(1, output.size());
		SoftClipEvidence ass = SoftClipEvidence.create(AES(ses), FWD, output.get(0));
		assertEquals(e.getBreakendSummary().localBreakend(), ass.getBreakendSummary());
		assertEquals("ACGTGGTCGACC", S(ass.getBreakendSequence()));
		assertEquals("", S(ass.getAnchorSequence()));
	}
	@Test
	public void should_call_simple_bwd_RP() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		MockSAMEvidenceSource ses = SES(100, 200);
		DiscordantReadPair e = (DiscordantReadPair)NRRP(ses, withSequence("ACGTGGTCGACC", DP(0, 450, "12M", false, 1, 1, "12M", true)));
		List<SAMRecord> output = go(pc, true, e);
		assertEquals(1, output.size());
		SoftClipEvidence ass = SoftClipEvidence.create(AES(ses), BWD, output.get(0));
		assertEquals(e.getBreakendSummary().localBreakend(), ass.getBreakendSummary());
		assertEquals("ACGTGGTCGACC", S(ass.getBreakendSequence()));
		assertEquals("", S(ass.getAnchorSequence()));
	}
	@Test
	public void should_handle_sc_repeated_kmers() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SoftClipEvidence sce = SCE(FWD, withSequence("TTTTTTTTTTTT", Read(0, 5, "6M6S")));
		List<SAMRecord> output = go(pc, true, sce);
		assertEquals(1, output.size());
		assertEquals("TTTTTTTTTTTT", S(output.get(0).getReadBases()));
	}
	@Test
	public void should_call_multiple_nonoverlapping_contigs() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SoftClipEvidence sce = SCE(FWD, withSequence("ACGTGGTCGACC", Read(0, 5, "6M6S")));
		DiscordantReadPair e = (DiscordantReadPair)NRRP(SES(100, 200), withSequence("GACCTCTACT", DP(0, 25, "10M", true, 1, 1, "10M", false)));
		List<SAMRecord> output = go(pc, true, sce, e);
		assertEquals(2, output.size());
		assertEquals("ACGTGGTCGACC", S(output.get(0).getReadBases()));
		assertEquals("NN" + "GACCTCTACT", S(output.get(1).getReadBases()));
	}
	@Test
	public void should_call_overlapping_sc_rp() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SoftClipEvidence sce = SCE(FWD, withSequence("ACGTGGTCGACC", Read(0, 50, "6M6S")));
		DiscordantReadPair e = (DiscordantReadPair)NRRP(SES(10, 200), withSequence("GACCTCTACT", DP(0, 25, "10M", true, 1, 1, "10M", false)));
		List<SAMRecord> output = go(pc, true, sce, e);
		assertEquals(1, output.size());
		assertEquals("ACGTGGTCGACCTCTACT", S(output.get(0).getReadBases()));
	}
	@Test
	public void should_remove_fully_reference_evidence() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SoftClipEvidence sce =     SCE(FWD, withSequence("ACGTGGTCGACC", Read(0, 50, "6M6S")));
		SoftClipEvidence fullref = SCE(FWD, withSequence("ACGTGG", Read(0, 50, "5M1S")));
		List<SAMRecord> output = go(pc, true, sce, fullref);
		assertEquals(1, output.size());
		assertEquals("ACGTGGTCGACC", S(output.get(0).getReadBases()));
		assertEquals(2, tracker.tracking_evidenceTotal());
		assertEquals(0, tracker.tracking_evidenceActive());
		assertEquals(0, caller.tracking_activeNodes());
	}
	@Test
	public void should_remove_fully_reference_evidence_before_end() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		int n = 8;
		DirectedEvidence[] input = new DirectedEvidence[2 * n];
		for (int i = 0; i < n; i++) {
			input[2*i] =      SCE(FWD, withSequence( "ACGTGGTCGACC", Read(0, 1000*i, "6M6S")));
			input[2*i + 1 ] = SCE(FWD, withSequence("AACGTGG", Read(0, 1000*i-1, "6M1S")));
		}
		NonReferenceContigAssembler caller = create(pc, true, input);
		for (int i = 0; i < n / 2; i++) {
			caller.next();
		}
		assertTrue(caller.tracking_activeNodes() <= 3);
	}
	@Test
	public void should_truncate_anchored_polyA_expansion_fwd() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SoftClipEvidence sce = SCE(FWD, withSequence("GTACAAAA", Read(0, 10, "4M4S")));
		DiscordantReadPair dp = (DiscordantReadPair)NRRP(SES(11, 100), withSequence("CAAAAT", DP(0, 1, "6M", true, 1, 1, "6M", false)));
		List<SAMRecord> output = go(pc, true, sce, dp);
		assertEquals(1, output.size());
		assertEquals("GTACAAAA", S(output.get(0).getReadBases()));
	}
	@Test
	public void should_truncate_anchored_polyA_expansion_bwd() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SoftClipEvidence sce = SCE(BWD, withSequence("AAAAGTAC", Read(0, 100, "4S4M")));
		DiscordantReadPair dp = (DiscordantReadPair)NRRP(SES(11, 100), withSequence("CAAAAG", DP(0, 110, "6M", false, 1, 1, "6M", true)));
		List<SAMRecord> output = go(pc, true, sce, dp);
		assertEquals(1, output.size());
		assertEquals("AAAAGTAC", S(output.get(0).getReadBases()));
	}
	@Test
	public void should_take_higher_weight_unanchored_assembly() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		pc.getAssemblyParameters().includePairAnchors = false;
		List<SAMRecord> output = go(pc, true,
				NRRP(SES(1, 100), withSequence("CAAAAAAAAAAAGT", DP(1, 100, "14M", false, 1, 1, "14M", true))),
				NRRP(SES(1, 100), withSequence("CAAAAAAAAAAGTC", DP(1, 101, "14M", false, 1, 1, "14M", true))));
		assertEquals(1, output.size());
		assertEquals("AAAAAAAAAAAGTC" + "NN", S(output.get(0).getReadBases()));
	}
	@Test
	public void should_handle_ambigious_bases() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		pc.getAssemblyParameters().includePairAnchors = false;
		List<SAMRecord> output = go(pc, true,
				NRRP(SES(1, 100), withSequence("AATAAGAGT", DP(1, 100, "9M", false, 1, 1, "9M", true))),
				NRRP(SES(1, 100), withSequence("AATGNGAGTC", DP(1, 101, "10M", false, 1, 1, "10M", true))));
		assertEquals(1, output.size());
		assertEquals("AATAAGAGTC" + "NN", S(output.get(0).getReadBases()));
	}
	@Test
	public void should_handle_multiple_candidate_offsets_for_a_single_kmer_position() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		pc.getAssemblyParameters().includePairAnchors = false;
		List<SAMRecord> output = go(pc, true,
				NRRP(SES(1, 100), withSequence("TAAAAA", DP(1, 100, "6M", false, 0, 1, "6M", true))));
		assertEquals(1, output.size());
		assertEquals("TAAAAA" + "NN", S(output.get(0).getReadBases()));
	}
	@Test
	public void fwd_anchor_length_calculation_should_not_include_any_nonreference_bases() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 2;
		List<SAMRecord> output = go(pc, true,
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 1, "9M1S"))),
				SCE(FWD, withSequence(        "ATTTTTTTTT", Read(0, 9, "1M9S"))));
		assertEquals(1, output.size());
		assertEquals("AAAAAAAAATTTTTTTTT", S(output.get(0).getReadBases()));
	}
	@Test
	public void bwd_anchor_length_calculation_should_not_include_any_nonreference_bases() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 2;
		List<SAMRecord> output = go(pc, true,
				SCE(BWD, withSequence("AAAAAAAAAT", Read(0, 10, "9S1M"))),
				SCE(BWD, withSequence(        "ATTTTTTTTT", Read(0, 10, "1S9M"))));
		assertEquals(1, output.size());
		assertEquals("AAAAAAAAATTTTTTTTT", S(output.get(0).getReadBases()));
	}
	@Test
	public void anchor_should_be_trimmed_at_max_length() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 3;
		List<SAMRecord> output = go(pc, true,
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 1, "9M1S"))),
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 2, "9M1S"))),
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 3, "9M1S"))),
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 4, "9M1S"))),
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 5, "9M1S"))),
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 6, "9M1S"))),
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 7, "9M1S"))),
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 8, "9M1S"))),
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 9, "9M1S"))),
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 10, "9M1S"))),
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 11, "9M1S"))),
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 12, "9M1S"))),
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 13, "9M1S"))),
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 14, "9M1S"))),
				SCE(FWD, withSequence(       "AATTTTTTTTT", Read(0, 18, "2M9S"))));
		assertEquals("AAAAAAAAAAATTTTTTTTT", S(output.get(0).getReadBases()));
	}
	@Test
	public void should_assembly_UM_into_anchor() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		pc.getAssemblyParameters().includePairAnchors = true;
		pc.getAssemblyParameters().pairAnchorMismatchIgnoreEndBases = 0;
		SAMRecord[] rp = OEA(0, 100, "10M", true);
		rp[0].setReadBases(B(S(RANDOM).substring(0, 10)));
		rp[1].setReadBases(B(S(RANDOM).substring(2, 12)));
		rp[1].setReadNegativeStrandFlag(true);
		NonReferenceReadPair nrrp = NRRP(SES(1, 100), rp);
		List<SAMRecord> output = go(pc, false, nrrp);
		assertEquals(S(RANDOM).substring(0, 12), S(output.get(0).getReadBases()));
		assertEquals(0, (int)output.get(0).getReferenceIndex());
		assertEquals(100, output.get(0).getAlignmentStart());
		assertEquals("10M2S", output.get(0).getCigarString());
	}
	@Test
	public void should_removeMisassembledPartialContigsDuringAssembly() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().removeMisassembledPartialContigsDuringAssembly = false;
		MockSAMEvidenceSource ses = new MockSAMEvidenceSource(getContext(), 50, 100);
		pc.getAssemblyParameters().k = 4;
		pc.getAssemblyParameters().maxExpectedBreakendLengthMultiple = 1.5f;
		List<DirectedEvidence> e = new ArrayList<>();
		for (int i = 1; i < 1000; i++) {
			e.add(NRRP(ses, OEA(0, i, "10M", true)));
		}
		List<SAMRecord> output = go(pc, false, e.toArray(new DirectedEvidence[0]));
		assertEquals(1, output.size());
		SAMRecord r = output.get(0);
		int sclenFull = SAMRecordUtil.getStartClipLength(r) + SAMRecordUtil.getEndClipLength(r);
		
		pc.getAssemblyParameters().removeMisassembledPartialContigsDuringAssembly = true;
		
		output = go(pc, false, e.toArray(new DirectedEvidence[0]));
		assertEquals(1, output.size());
		r = output.get(0);
		int sclenTruncated = SAMRecordUtil.getStartClipLength(r) + SAMRecordUtil.getEndClipLength(r);
		Assert.assertNotEquals(sclenFull, sclenTruncated);
	}
	@Test
	public void should_call_suboptimal_contigs_if_graph_size_limit_exceeded() {
		ProcessingContext pc = getContext();
		MockSAMEvidenceSource ses = SES(10, 10);
		pc.getAssemblyParameters().k = 4;
		pc.getAssemblyParameters().maxExpectedBreakendLengthMultiple = 1;
		pc.getAssemblyParameters().positional.retainWidthMultiple = 1;
		pc.getAssemblyParameters().positional.flushWidthMultiple = 2;
		List<DirectedEvidence> e = new ArrayList<>();
		for (int i = 0; i < 100; i++) {
			for (int j = 0; j < i + 1; j++) {
				e.add(SCE(FWD, ses, withReadName(String.format("%d-%d", i, j), withSequence("AACCGGTT", Read(0, i, "4M4S")))[0]));
			}
		}
		List<SAMRecord> output = go(pc, false, e.toArray(new DirectedEvidence[0]));
		assertEquals(100, output.size());
		// unbounded execution would call high->low thus call the last contig first
		// bounded execution will force the early calls to be made first
		List<Integer> startPos = output.stream().map(x -> x.getAlignmentStart()).collect(Collectors.toList());
		assertFalse(Ordering.natural().reverse().isOrdered(startPos));
	}
	@Test
	public void should_not_assembly_when_graph_density_exceeds_maximum() {
		ProcessingContext pc = getContext();
		MockSAMEvidenceSource ses = SES(10, 10);
		pc.getAssemblyParameters().k = 4;
		pc.getAssemblyParameters().maxExpectedBreakendLengthMultiple = 1;
		pc.getAssemblyParameters().positional.maximumNodeDensity = 0.1f;
		List<DirectedEvidence> e = new ArrayList<>();
		for (int i = 0; i < 100; i++) {
			e.add(SCE(FWD, ses, withReadName(String.format("%d-%d", i, 0), withSequence("AAAATTGG", Read(0, i, "4M4S")))[0]));
			e.add(SCE(FWD, ses, withReadName(String.format("%d-%d", i, 1), withSequence("AAAACCGG", Read(0, i, "4M4S")))[0]));
		}
		List<SAMRecord> output = go(pc, false, e.toArray(new DirectedEvidence[0]));
		// the final advance doesn't get flushed since the underlying stream
		// is considered to advance to Integer.MAX_VALUE at end of input
		assertTrue(output.size() <= 20);
		
		// shouldn't get throttled
		pc.getAssemblyParameters().positional.maximumNodeDensity = 10;
		output = go(pc, false, e.toArray(new DirectedEvidence[0]));
		assertEquals(2 * 100, output.size());
	}
	@Test
	public void should_remove_misassembled_partial_paths() {
		ProcessingContext pc = getContext();
		MockSAMEvidenceSource ses = SES(10, 10);
		pc.getAssemblyParameters().k = 4;
		pc.getAssemblyParameters().maxExpectedBreakendLengthMultiple = 1;
		List<DirectedEvidence> e = new ArrayList<>();
		for (int i = 1; i < 100; i++) {
			e.add(SCE(FWD, ses, withReadName(String.format("%d-%d", 0, i), withSequence("TAAAAAAAAAAAAAAAAAAAAAAAAAAAA", Read(0, i*20, "4M25S")))[0]));
		}
		List<SAMRecord> output = go(pc, false, e.toArray(new DirectedEvidence[0]));
		// the rest should have been removed
		assertEquals(1, output.size());
	}
	@Test
	public void should_not_remove_misassembled_partial_paths_when_could_be_collapsed_with_adjacent_read() {
		ProcessingContext pc = getContext();
		MockSAMEvidenceSource ses = SES(10, 10);
		pc.getAssemblyParameters().k = 2;
		pc.getAssemblyParameters().maxExpectedBreakendLengthMultiple = 5;
		pc.getAssemblyParameters().removeMisassembledPartialContigsDuringAssembly = true;
		pc.getAssemblyParameters().positional.flushWidthMultiple = 100000;
		pc.getAssemblyParameters().positional.maxPathLengthMultiple = 4;
		List<DirectedEvidence> e = new ArrayList<>();
		for (int i = 1; i < 1000; i += 1) {
			SAMRecord[] rp = OEA(0, i, "2M", true);
			rp[0].setReadBases(B("GG"));
			rp[1].setReadBases(B("AA")); // expect reverse comp of read sequence
			e.add(NRRP(ses, rp));
		}
		go(pc, true, e.toArray(new DirectedEvidence[0]));
	}
}
