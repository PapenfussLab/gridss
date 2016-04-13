package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

import org.junit.Test;

import au.edu.wehi.idsv.AssemblyEvidence;
import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.DiscordantReadPair;
import au.edu.wehi.idsv.NonReferenceReadPair;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.SAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.SoftClipEvidence;
import au.edu.wehi.idsv.SpanningSAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.TestHelper;

import com.google.common.collect.Lists;


public class NonReferenceContigAssemblerTest extends TestHelper {
	private NonReferenceContigAssembler caller;
	private EvidenceTracker tracker;
	public List<SAMRecordAssemblyEvidence> go(ProcessingContext pc, boolean collapse, DirectedEvidence... input) {
		caller = create(pc, collapse, input);
		List<SAMRecordAssemblyEvidence> assemblies = Lists.newArrayList(caller);
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
		SupportNodeIterator supportIt = new SupportNodeIterator(k, Arrays.stream(input).iterator(), aes.getMaxConcordantFragmentSize(), tracker, pc.getAssemblyParameters().includePairAnchors);
		AggregateNodeIterator agIt = new AggregateNodeIterator(supportIt);
		Iterator<KmerPathNode> pnIt = new PathNodeIterator(agIt, maxPathLength, k);
		if (collapse) {
			pnIt = new PathCollapseIterator(pnIt, k, maxPathCollapseLength, pc.getAssemblyParameters().errorCorrection.maxBaseMismatchForCollapse, pc.getAssemblyParameters().errorCorrection.collapseBubblesOnly, 0);
			pnIt = new PathSimplificationIterator(pnIt, maxPathLength, maxEvidenceWidth);
		}
		caller = new NonReferenceContigAssembler(pnIt, 0, maxEvidenceWidth + maxReadLength + 2, maxReadLength, k, aes, tracker, "test");
		return caller;
	}
	@Test
	public void should_call_simple_fwd_SC() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SoftClipEvidence sce = SCE(FWD, withSequence("ACGTGGTCGACC", Read(0, 5, "6M6S")));
		List<SAMRecordAssemblyEvidence> output = go(pc, true, sce);
		assertEquals(1, output.size());
		assertEquals("ACGTGGTCGACC", S(output.get(0).getAssemblySequence()));
		assertEquals(new BreakendSummary(0, FWD, 10, 10), output.get(0).getBreakendSummary());
		assertEquals(6, output.get(0).getAssemblyAnchorLength());
	}
	@Test
	public void should_call_simple_bwd_SC() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SoftClipEvidence sce = SCE(BWD, withSequence("ACGTGGTCGACC", Read(0, 10, "6S6M")));
		List<SAMRecordAssemblyEvidence> output = go(pc, true, sce);
		assertEquals(1, output.size());
		assertEquals("ACGTGGTCGACC", S(output.get(0).getAssemblySequence()));
		assertEquals(new BreakendSummary(0, BWD, 10, 10), output.get(0).getBreakendSummary());
		assertEquals(6, output.get(0).getAssemblyAnchorLength());
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
		List<SAMRecordAssemblyEvidence> output = go(pc, true, fwd, bwd);
		assertEquals(1, output.size());
		AssemblyEvidence e = output.get(0).getSpannedIndels().get(0);
		assertEquals("ACGTGGCTCGACC", S(e.getAssemblySequence()));
		assertEquals(new BreakpointSummary(0, FWD, 6, 6, 0, BWD, 8, 8), e.getBreakendSummary());
		assertEquals("C", ((SpanningSAMRecordAssemblyEvidence)e).getUntemplatedSequence());
	}
	@Test
	public void should_call_simple_fwd_RP() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		DiscordantReadPair e = (DiscordantReadPair)NRRP(SES(100, 200), withSequence("ACGTGGTCGACC", DP(0, 50, "12M", true, 1, 1, "12M", false)));
		List<SAMRecordAssemblyEvidence> output = go(pc, true, e);
		assertEquals(1, output.size());
		assertEquals("ACGTGGTCGACC", S(output.get(0).getAssemblySequence()));
		assertEquals(e.getBreakendSummary().localBreakend(), output.get(0).getBreakendSummary());
	}
	@Test
	public void should_call_simple_bwd_RP() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		DiscordantReadPair e = (DiscordantReadPair)NRRP(SES(100, 200), withSequence("ACGTGGTCGACC", DP(0, 450, "12M", false, 1, 1, "12M", true)));
		List<SAMRecordAssemblyEvidence> output = go(pc, true, e);
		assertEquals(1, output.size());
		assertEquals("ACGTGGTCGACC", S(output.get(0).getAssemblySequence()));
		assertEquals(e.getBreakendSummary().localBreakend(), output.get(0).getBreakendSummary());
	}
	@Test
	public void should_handle_sc_repeated_kmers() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SoftClipEvidence sce = SCE(FWD, withSequence("TTTTTTTTTTTT", Read(0, 5, "6M6S")));
		List<SAMRecordAssemblyEvidence> output = go(pc, true, sce);
		assertEquals(1, output.size());
		assertEquals("TTTTTTTTTTTT", S(output.get(0).getAssemblySequence()));
	}
	@Test
	public void should_call_multiple_nonoverlapping_contigs() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SoftClipEvidence sce = SCE(FWD, withSequence("ACGTGGTCGACC", Read(0, 5, "6M6S")));
		DiscordantReadPair e = (DiscordantReadPair)NRRP(SES(100, 200), withSequence("GACCTCTACT", DP(0, 25, "10M", true, 1, 1, "10M", false)));
		List<SAMRecordAssemblyEvidence> output = go(pc, true, sce, e);
		assertEquals(2, output.size());
		assertEquals("ACGTGGTCGACC", S(output.get(0).getAssemblySequence()));
		assertEquals("GACCTCTACT", S(output.get(1).getAssemblySequence()));
		assertEquals(sce.getBreakendSummary(), output.get(0).getBreakendSummary());
		assertEquals(e.getBreakendSummary().localBreakend(), output.get(1).getBreakendSummary());
	}
	@Test
	public void should_call_overlapping_sc_rp() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SoftClipEvidence sce = SCE(FWD, withSequence("ACGTGGTCGACC", Read(0, 50, "6M6S")));
		DiscordantReadPair e = (DiscordantReadPair)NRRP(SES(10, 200), withSequence("GACCTCTACT", DP(0, 25, "10M", true, 1, 1, "10M", false)));
		List<SAMRecordAssemblyEvidence> output = go(pc, true, sce, e);
		assertEquals(1, output.size());
		assertEquals("ACGTGGTCGACCTCTACT", S(output.get(0).getAssemblySequence()));
	}
	@Test
	public void should_remove_fully_reference_evidence() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SoftClipEvidence sce =     SCE(FWD, withSequence("ACGTGGTCGACC", Read(0, 50, "6M6S")));
		SoftClipEvidence fullref = SCE(FWD, withSequence("ACGTGG", Read(0, 50, "5M1S")));
		List<SAMRecordAssemblyEvidence> output = go(pc, true, sce, fullref);
		assertEquals(1, output.size());
		assertEquals("ACGTGGTCGACC", S(output.get(0).getAssemblySequence()));
		assertEquals(2, tracker.tracking_evidenceTotal());
		assertEquals(0, tracker.tracking_evidenceActive());
		assertEquals(0, caller.tracking_activeNodes());
	}
	@Test
	public void should_remove_fully_reference_evidence_before_end() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		int n = 64;
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
		List<SAMRecordAssemblyEvidence> output = go(pc, true, sce, dp);
		assertEquals(1, output.size());
		assertEquals("GTACAAAA", S(output.get(0).getAssemblySequence()));
	}
	@Test
	public void should_truncate_anchored_polyA_expansion_bwd() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SoftClipEvidence sce = SCE(BWD, withSequence("AAAAGTAC", Read(0, 100, "4S4M")));
		DiscordantReadPair dp = (DiscordantReadPair)NRRP(SES(11, 100), withSequence("CAAAAG", DP(0, 110, "6M", false, 1, 1, "6M", true)));
		List<SAMRecordAssemblyEvidence> output = go(pc, true, sce, dp);
		assertEquals(1, output.size());
		assertEquals("AAAAGTAC", S(output.get(0).getAssemblySequence()));
	}
	@Test
	public void should_take_higher_weight_unanchored_assembly() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		pc.getAssemblyParameters().includePairAnchors = false;
		List<SAMRecordAssemblyEvidence> output = go(pc, true,
				NRRP(SES(1, 100), withSequence("CAAAAAAAAAAAGT", DP(1, 100, "14M", false, 1, 1, "14M", true))),
				NRRP(SES(1, 100), withSequence("CAAAAAAAAAAGTC", DP(1, 101, "14M", false, 1, 1, "14M", true))));
		assertEquals(1, output.size());
		assertEquals("AAAAAAAAAAAGTC", S(output.get(0).getAssemblySequence()));
	}
	@Test
	public void should_handle_ambigious_bases() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		pc.getAssemblyParameters().includePairAnchors = false;
		List<SAMRecordAssemblyEvidence> output = go(pc, true,
				NRRP(SES(1, 100), withSequence("CAAAAAAAAAAAGT", DP(1, 100, "14M", false, 1, 1, "14M", true))),
				NRRP(SES(1, 100), withSequence("CAAAAANAAAAGTC", DP(1, 101, "14M", false, 1, 1, "14M", true))));
		assertEquals(1, output.size());
		assertEquals("AAAAAAAAAAAGTC", S(output.get(0).getAssemblySequence()));
	}
	@Test
	public void should_handle_multiple_candidate_offsets_for_a_single_kmer_position() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		pc.getAssemblyParameters().includePairAnchors = false;
		List<SAMRecordAssemblyEvidence> output = go(pc, true,
				NRRP(SES(1, 100), withSequence("TAAAAA", DP(1, 100, "6M", false, 0, 1, "6M", true))));
		assertEquals(1, output.size());
		assertEquals("TAAAAA", S(output.get(0).getAssemblySequence()));
	}
	@Test
	public void fwd_anchor_length_calculation_should_not_include_any_nonreference_bases() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 2;
		List<SAMRecordAssemblyEvidence> output = go(pc, true,
				SCE(FWD, withSequence("AAAAAAAAAT", Read(0, 1, "9M1S"))),
				SCE(FWD, withSequence(        "ATTTTTTTTT", Read(0, 9, "1M9S"))));
		assertEquals(1, output.size());
		assertEquals("AAAAAAAAATTTTTTTTT", S(output.get(0).getAssemblySequence()));
	}
	@Test
	public void bwd_anchor_length_calculation_should_not_include_any_nonreference_bases() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 2;
		List<SAMRecordAssemblyEvidence> output = go(pc, true,
				SCE(BWD, withSequence("AAAAAAAAAT", Read(0, 10, "9S1M"))),
				SCE(BWD, withSequence(        "ATTTTTTTTT", Read(0, 10, "1S9M"))));
		assertEquals(1, output.size());
		assertEquals("AAAAAAAAATTTTTTTTT", S(output.get(0).getAssemblySequence()));
	}
	@Test
	public void anchor_should_be_trimmed_at_max_length() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 3;
		List<SAMRecordAssemblyEvidence> output = go(pc, true,
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
		assertEquals("AAAAAAAAAAATTTTTTTTT", S(output.get(0).getAssemblySequence()));
	}
	@Test
	public void should_assembly_UM_into_anchor() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		SAMRecord[] rp = OEA(0, 100, "10M", true);
		rp[0].setReadBases(B(S(RANDOM).substring(0, 10)));
		rp[1].setReadBases(B(S(RANDOM).substring(2, 12)));
		rp[1].setReadNegativeStrandFlag(true);
		NonReferenceReadPair nrrp = NRRP(SES(1, 100), rp);
		List<SAMRecordAssemblyEvidence> output = go(pc, false, nrrp);
		assertEquals(S(RANDOM).substring(0, 12), S(output.get(0).getAssemblySequence()));
		assertEquals(new BreakendSummary(0, FWD, 100 + 10 - 1, 100 + 10 - 1), output.get(0).getBreakendSummary());
	}
	@Test
	public void should_strip_evidence_forming_contig_longer_than_maxExpectedBreakendLengthMultiple() {
		ProcessingContext pc = getContext();
		MockSAMEvidenceSource ses = SES(50, 100);
		pc.getAssemblyParameters().k = 4;
		List<DirectedEvidence> e = new ArrayList<>();
		e.add(SCE(FWD, withSequence("ACGTAACCGGTT", Read(0, 1, "6M6S"))));
		for (int i = 1; i < 1000; i++) {
			e.add(NRRP(ses, OEA(0, i, "10M", true)));
		}
		List<SAMRecordAssemblyEvidence> output = go(pc, false, e.toArray(new DirectedEvidence[0]));
		assertEquals(1, output.size());
	}
}
