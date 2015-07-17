package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

import org.junit.Test;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.DiscordantReadPair;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.SAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.SmallIndelSAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.SoftClipEvidence;
import au.edu.wehi.idsv.TestHelper;

import com.google.common.collect.Lists;


public class NonReferenceContigAssemblerTest extends TestHelper {
	private NonReferenceContigAssembler caller;
	private EvidenceTracker trackedIt;
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
		int maxPathLength = pc.getAssemblyParameters().positionalMaxPathLengthInBases(maxReadLength);
		int maxPathCollapseLength = pc.getAssemblyParameters().positionalMaxPathCollapseLengthInBases(maxReadLength);
		SupportNodeIterator supportIt = new SupportNodeIterator(k, Arrays.stream(input).iterator(), aes.getMaxConcordantFragmentSize());
		trackedIt = new EvidenceTracker(supportIt);
		AggregateNodeIterator agIt = new AggregateNodeIterator(trackedIt);
		Iterator<KmerPathNode> pnIt = new PathNodeIterator(agIt, maxPathLength, k);
		if (collapse) {
			pnIt = new PathCollapseIterator(pnIt, k, maxPathCollapseLength, pc.getAssemblyParameters().maxBaseMismatchForCollapse, pc.getAssemblyParameters().collapseBubblesOnly, 0);
			pnIt = new PathSimplificationIterator(pnIt, maxPathLength, maxEvidenceWidth);
		}
		caller = new NonReferenceContigAssembler(pnIt, 0, maxEvidenceWidth, maxReadLength, k, aes, trackedIt);
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
		assertTrue(output.get(0) instanceof SmallIndelSAMRecordAssemblyEvidence);
		assertEquals("ACGTGGCTCGACC", S(output.get(0).getAssemblySequence()));
		assertEquals(new BreakpointSummary(0, FWD, 6, 6, 0, BWD, 8, 8), output.get(0).getBreakendSummary());
		assertEquals("C", ((SmallIndelSAMRecordAssemblyEvidence)output.get(0)).getUntemplatedSequence());
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
		SoftClipEvidence sce = SCE(FWD, withSequence("ACGTGGTCGACC", Read(0, 50, "6M6S")));
		SoftClipEvidence fullref = SCE(FWD, withSequence("ACGTGG", Read(0, 50, "5M1S")));
		List<SAMRecordAssemblyEvidence> output = go(pc, true, sce, fullref);
		assertEquals(1, output.size());
		assertEquals("ACGTGGTCGACC", S(output.get(0).getAssemblySequence()));
		assertEquals(0, trackedIt.tracking_supportNodeCount());
		assertEquals(0, caller.tracking_activeNodes());
	}
	@Test
	public void should_remove_fully_reference_evidence_before_end() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().k = 4;
		int n = 16;
		DirectedEvidence[] input = new DirectedEvidence[2 * n];
		for (int i = 0; i < n; i++) {
			input[2*i] = SCE(FWD, withSequence("ACGTGGTCGACC", Read(0, 1000*i, "6M6S")));
			input[2*i + 1 ] = SCE(FWD, withSequence("AACGTGG", Read(0, 1000*i-1, "6M1S")));
		}
		NonReferenceContigAssembler caller = create(pc, true, input);
		for (int i = 0; i < n / 2; i++) {
			caller.next();
		}
		assertTrue(caller.tracking_activeNodes() <= 3);
	}
}
