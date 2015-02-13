package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;



public class SequentialEvidenceAnnotatorTest extends TestHelper {
	public VariantContextDirectedEvidence go(List<DirectedEvidence> evidence, VariantContextDirectedEvidence toAnnotate) {
		return Lists.newArrayList(new SequentialEvidenceAnnotator(getContext(), ImmutableList.of(toAnnotate).iterator(), evidence.iterator(), 300, true, null)).get(0);
	}
	@Test
	public void should_add_breakpoint_supporting_evidence() {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext());
		builder
			.phredScore(1)
			.loc("polyA", 1, 1)
			.alleles("A", "A[polyA:10[");
		VariantContextDirectedBreakpoint result = (VariantContextDirectedBreakpoint)go(L(
				(DirectedEvidence)SoftClipEvidence.create(SES(), FWD, Read(0, 1, "1M3S"), Read(0, 10, "3M"))
			), (VariantContextDirectedEvidence)builder.make());
		assertEquals(1, result.getBreakpointEvidenceCountSoftClip(null));
	}
	@Test
	public void should_add_breakend_supporting_evidence() {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext());
		builder
			.phredScore(1)
			.loc("polyA", 1, 1)
			.alleles("A", "A[polyA:10[");
		VariantContextDirectedEvidence result = go(L(
				(DirectedEvidence)SoftClipEvidence.create(SES(), FWD, Read(0, 1, "1M3S"))
			), (VariantContextDirectedEvidence)builder.make());
		assertEquals(1, result.getBreakendEvidenceCountSoftClip(null));
	}
	@Test
	public void should_use_breakpoint_assembly_sequence() {
		SAMRecordAssemblyEvidence be = AssemblyFactory.createAnchored(getContext(), AES(), FWD, Sets.<DirectedEvidence>newHashSet(),
				0, 1, 1, B("TGGAAT"), B("TAAAAT"), 0, 0);
		// Untemplated sequence
		SAMRecord r = Read(0, 10, "2S3M");
		r.setReadBases(B("GGAAT"));
		RealignedSAMRecordAssemblyEvidence ass = (RealignedSAMRecordAssemblyEvidence) AssemblyFactory.incorporateRealignment(getContext(), be, r);
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext());
		builder
			.phredScore(1)
			.loc("polyA", 0, 1)
			.alleles("A", "A[polyA:10[");
		VariantContextDirectedBreakpoint result = (VariantContextDirectedBreakpoint)go(L(
				(DirectedEvidence)ass
			), (VariantContextDirectedEvidence)builder.make());
		assertEquals("AGG[polyA:10[", result.getAlternateAllele(0).getDisplayString());
	}
	@Test
	public void should_merge_supporting_evidence() {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext());
		builder
			.phredScore(1)
			.loc("polyA", 1, 1)
			.alleles("A", "A[polyA:10[");
		VariantContextDirectedBreakpoint result = (VariantContextDirectedBreakpoint)go(L(
				(DirectedEvidence)SoftClipEvidence.create(SES(), FWD, Read(0, 1, "1M2S"), Read(0, 10, "3M")),
				(DirectedEvidence)SoftClipEvidence.create(SES(), FWD, Read(0, 1, "1M3S"), Read(0, 10, "3M"))
			), (VariantContextDirectedEvidence)builder.make());
		assertEquals(2, result.getBreakpointEvidenceCountSoftClip(null));
	}
	@Test
	public void should_incorporate_sc_margin_when_adding_evidence() {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext());
		builder
			.phredScore(1)
			.loc("polyA", 2, 2)
			.alleles("A", "A[polyA:10[");
		VariantContextDirectedBreakpoint result = (VariantContextDirectedBreakpoint)go(L(
				(DirectedEvidence)SoftClipEvidence.create(SES(), FWD, Read(0, 1, "1M3S"))
			), (VariantContextDirectedEvidence)builder.make());
		assertEquals(1, result.getBreakendEvidenceCountSoftClip(null));
	}
	@Test
	public void should_output_variants_in_order() {
		List<VariantContextDirectedEvidence> calls = new ArrayList<VariantContextDirectedEvidence>();
		calls.add((VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
				breakpoint(new BreakpointSummary(0, FWD, 1, 3, 0, BWD, 10, 10), "GTAC");
				phredScore(10);
			}}.make());
		calls.add((VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
				breakpoint(new BreakpointSummary(0, FWD, 2, 4, 0, BWD, 10, 10), "GTAC");
				phredScore(20);
			}}.make());
		List<DirectedEvidence> evidence = new ArrayList<DirectedEvidence>();
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(new SequentialEvidenceAnnotator(getContext(), calls.iterator(), evidence.iterator(), 300, true, null));
		assertEquals(1, result.get(0).getBreakendSummary().start);
		assertEquals(2, result.get(1).getBreakendSummary().start);
	}
	@Test
	public void should_assign_supporting_evidence_to_a_single_variant() {
		List<VariantContextDirectedEvidence> calls = new ArrayList<VariantContextDirectedEvidence>();
		calls.add((VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
				breakpoint(new BreakpointSummary(0, FWD, 11, 12, 0, BWD, 10, 10), "");
				phredScore(10);
			}}.make());
		calls.add((VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
				breakpoint(new BreakpointSummary(0, FWD, 12, 14, 0, BWD, 10, 10), "");
				phredScore(20);
			}}.make());
		List<DirectedEvidence> evidence = new ArrayList<DirectedEvidence>();
		evidence.add(SoftClipEvidence.create(SES(), FWD, withSequence("NNNNN", Read(0, 12, "1M4S"))[0], withSequence("NNN", Read(0, 10, "3M"))[0]));
		evidence.add(SoftClipEvidence.create(SES(), FWD, withSequence("NNNN", Read(0, 12, "1M3S"))[0], withSequence("NNN", Read(0, 10, "3M"))[0]));
		
		BreakpointSummary bp = new BreakpointSummary(0,  FWD, 12, 12, 0, BWD, 10, 10); 
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(new SequentialEvidenceAnnotator(getContext(), calls.iterator(), evidence.iterator(), 300, true, null));
		assertNotEquals(bp, result.get(0).getBreakendSummary());
		assertEquals(bp, result.get(1).getBreakendSummary());
		assertEquals(0, ((VariantContextDirectedBreakpoint)result.get(0)).getBreakpointEvidenceCountSoftClip(null));
		assertEquals(2, ((VariantContextDirectedBreakpoint)result.get(1)).getBreakpointEvidenceCountSoftClip(null));
		
		result = Lists.newArrayList(new SequentialEvidenceAnnotator(getContext(), calls.iterator(), evidence.iterator(), 300, false, null));
		assertEquals(bp, result.get(0).getBreakendSummary());
		assertEquals(bp, result.get(1).getBreakendSummary());
		assertEquals(2, ((VariantContextDirectedBreakpoint)result.get(0)).getBreakpointEvidenceCountSoftClip(null));
		assertEquals(2, ((VariantContextDirectedBreakpoint)result.get(1)).getBreakpointEvidenceCountSoftClip(null));
	}
	@Test
	public void should_not_add_breakend_nonsupporting_evidence() {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext());
		builder
			.phredScore(1)
			.loc("polyA", 1, 1)
			.alleles("A", "A[polyA:10[");
		VariantContextDirectedBreakpoint result = (VariantContextDirectedBreakpoint)go(L(
				(DirectedEvidence)SoftClipEvidence.create(SES(), BWD, Read(0, 1, "3S1M"), Read(0, 10, "3M"))
			), (VariantContextDirectedEvidence)builder.make());
		assertEquals(0, result.getBreakpointEvidenceCountSoftClip(null));
		assertEquals(0, result.getBreakendEvidenceCountSoftClip(null));
	}
	@Test
	public void should_not_add_breakpoint_nonsupporting_evidence() {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext());
		builder
			.phredScore(1)
			.loc("polyA", 1, 1)
			.alleles("A", "A[polyA:10[");
		VariantContextDirectedBreakpoint result = (VariantContextDirectedBreakpoint)go(L(
				(DirectedEvidence)SoftClipEvidence.create(SES(), FWD, withSequence("TTTTTTTT", Read(0, 5, "3S1M3S"))[0], withSequence("TTT", Read(0, 22, "3M"))[0])
			), (VariantContextDirectedEvidence)builder.make());
		assertEquals(0, result.getBreakpointEvidenceCountSoftClip(null));
		assertEquals(0, result.getBreakendEvidenceCountSoftClip(null));
	}
	@Test
	public void should_allocate_in_order() {
		int fragSize = 8;
		List<DirectedEvidence> evidence = new ArrayList<DirectedEvidence>();
		List<VariantContextDirectedBreakpoint> calls = new ArrayList<VariantContextDirectedBreakpoint>();
		ProcessingContext pc = getContext();
		MockSAMEvidenceSource ses = SES(fragSize, fragSize);
		int testSize = 16;
		for (int i = 1; i <= testSize; i++) {
			for (int j = 1; j <= testSize; j++) {
				IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(pc);
				builder
					.phredScore(1)
					.breakpoint(new BreakpointSummary(0, FWD, i, i, 1, BWD, j, j), "")
					.attribute("EVENT", String.format("%d-%d", i, j));
				calls.add((VariantContextDirectedBreakpoint)builder.make());
				// And the matching reverse
				IdsvVariantContextBuilder builder2 = new IdsvVariantContextBuilder(pc);
				builder2
					.phredScore(1)
					.breakpoint(new BreakpointSummary(1, BWD, j, j, 0, FWD, i, i), "")
					.attribute("EVENT", String.format("%d-%d", i, j));
				calls.add((VariantContextDirectedBreakpoint)builder2.make());
			}
		}
		for (int i = 1; i <= testSize; i++) {
			for (int j = 1; j <= testSize; j++) {
				SAMRecord[] dp = DP(0, i, "1M", true, 1, j, "1M", false);
				evidence.add(NonReferenceReadPair.create(dp[0], dp[1], ses));
				evidence.add(NonReferenceReadPair.create(dp[1], dp[0], ses));
			}
		}
		Collections.sort(calls, DirectedEvidenceOrder.ByNatural);
		Collections.sort(evidence, DirectedEvidenceOrder.ByNatural);
		ArrayList<VariantContextDirectedEvidence> result = Lists.newArrayList(new SequentialEvidenceAnnotator(pc, calls.iterator(), evidence.iterator(), fragSize, true, null));
		assertEquals(calls.size(), result.size());
		
		double expectedEvidence = 0;
		for (DirectedEvidence e : evidence) {
			DiscordantReadPair bp = (DiscordantReadPair) e;
			expectedEvidence += bp.getBreakpointQual();
		}
		double annotatedEvidence = 0;
		for (IdsvVariantContext e : result) {
			annotatedEvidence += e.getPhredScaledQual();
		}
		// each piece of evidence should be assigned to a single breakpoint so totals should match
		assertEquals(expectedEvidence, annotatedEvidence, 0.1);
	}
	@Test
	public void should_allocate_related_evidence_to_different_breakends_when_breakpoint_self_overlapping() {
		List<VariantContextDirectedBreakpoint> calls = new ArrayList<VariantContextDirectedBreakpoint>();
		List<DirectedEvidence> evidence = new ArrayList<DirectedEvidence>();
		ProcessingContext pc = getContext();
		MockSAMEvidenceSource ses = SES(10, 10);
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(pc);
		builder
			.phredScore(1)
			.breakpoint(new BreakpointSummary(0, FWD, 1, 10, 0, FWD, 1, 10), "")
			.attribute("EVENT", "SELFINTERSECTING")
			.id("low")
			.attribute("MATEID", "high");
		calls.add((VariantContextDirectedBreakpoint)builder.make());
		builder = new IdsvVariantContextBuilder(pc);
		builder
			.phredScore(1)
			.breakpoint(new BreakpointSummary(0, FWD, 2, 10, 0, FWD, 2, 10), "")
			.attribute("EVENT", "SELFINTERSECTING")
			.id("high")
			.attribute("MATEID", "low");
		calls.add((VariantContextDirectedBreakpoint)builder.make());
		
		SAMRecord[] dp = DP(0, 5, "1M", true, 0, 4, "1M", true);
		evidence.add(NonReferenceReadPair.create(dp[0], dp[1], ses));
		evidence.add(NonReferenceReadPair.create(dp[1], dp[0], ses));
		
		Collections.sort(calls, DirectedEvidenceOrder.ByNatural);
		Collections.sort(evidence, DirectedEvidenceOrder.ByNatural);
		ArrayList<VariantContextDirectedEvidence> result = Lists.newArrayList(new SequentialEvidenceAnnotator(pc, calls.iterator(), evidence.iterator(), 10, true, null));
		assertEquals(calls.size(), result.size());
		assertEquals(1, ((VariantContextDirectedBreakpoint)result.get(0)).getBreakpointEvidenceCountReadPair(EvidenceSubset.ALL));
		assertEquals(1, ((VariantContextDirectedBreakpoint)result.get(1)).getBreakpointEvidenceCountReadPair(EvidenceSubset.ALL));
	}
	private List<VariantContextDirectedBreakpoint> buildPair(ProcessingContext pc, BreakpointSummary bs, String eventid) {
		return buildPair(pc, bs, eventid, 1);
	}
	private List<VariantContextDirectedBreakpoint> buildPair(ProcessingContext pc, BreakpointSummary bs, String eventid, float score) {
		List<VariantContextDirectedBreakpoint> list = new ArrayList<VariantContextDirectedBreakpoint>();
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(pc);
		builder
			.phredScore(score)
			.breakpoint(bs, "")
			.attribute("EVENT", eventid)
			.id("low" + eventid)
			.attribute("MATEID", "high" + eventid);
		list.add((VariantContextDirectedBreakpoint)builder.make());
		builder = new IdsvVariantContextBuilder(pc);
		builder
			.phredScore(score)
			.breakpoint(bs.remoteBreakpoint(), "")
			.attribute("EVENT", eventid)
			.id("high" + eventid)
			.attribute("MATEID", "low" + eventid);
		list.add((VariantContextDirectedBreakpoint)builder.make());
		return list;
	}
	@Test
	public void should_disambiguate_equally_weighted_calls_by_low_high_breakend_coordinates() {
		List<VariantContextDirectedBreakpoint> calls = new ArrayList<VariantContextDirectedBreakpoint>();
		List<DirectedEvidence> evidence = new ArrayList<DirectedEvidence>();
		ProcessingContext pc = getContext();
		MockSAMEvidenceSource ses = SES(10, 10);
		
		calls.addAll(buildPair(pc, new BreakpointSummary(0, FWD, 1, 10, 1, FWD, 11, 20), "1"));
		calls.addAll(buildPair(pc, new BreakpointSummary(0, FWD, 11, 20, 1, FWD, 1, 10), "2"));
		
		SAMRecord[] dp = DP(0, 1, "1M", true, 1, 1, "1M", true);
		evidence.add(NonReferenceReadPair.create(dp[0], dp[1], ses));
		evidence.add(NonReferenceReadPair.create(dp[1], dp[0], ses));
		
		Collections.sort(calls, DirectedEvidenceOrder.ByNatural);
		Collections.sort(evidence, DirectedEvidenceOrder.ByNatural);
		ArrayList<VariantContextDirectedEvidence> result = Lists.newArrayList(new SequentialEvidenceAnnotator(pc, calls.iterator(), evidence.iterator(), 20, true, null));
		assertEquals(calls.size(), result.size());
		assertEquals(1, ((VariantContextDirectedBreakpoint)result.get(0)).getBreakpointEvidenceCountReadPair(EvidenceSubset.ALL)); // low1
		assertEquals(0, ((VariantContextDirectedBreakpoint)result.get(1)).getBreakpointEvidenceCountReadPair(EvidenceSubset.ALL)); // low2
		assertEquals(0, ((VariantContextDirectedBreakpoint)result.get(2)).getBreakpointEvidenceCountReadPair(EvidenceSubset.ALL)); // high2
		assertEquals(1, ((VariantContextDirectedBreakpoint)result.get(3)).getBreakpointEvidenceCountReadPair(EvidenceSubset.ALL)); // high1
	}
	@Test
	public void should_disambiguate_calls_by_weight() {
		List<VariantContextDirectedBreakpoint> calls = new ArrayList<VariantContextDirectedBreakpoint>();
		List<DirectedEvidence> evidence = new ArrayList<DirectedEvidence>();
		ProcessingContext pc = getContext();
		MockSAMEvidenceSource ses = SES(10, 10);
		
		calls.addAll(buildPair(pc, new BreakpointSummary(0, FWD, 1, 10, 1, FWD, 11, 20), "1", 1));
		calls.addAll(buildPair(pc, new BreakpointSummary(0, FWD, 11, 20, 1, FWD, 1, 10), "2", 2));
		
		SAMRecord[] dp = DP(0, 1, "1M", true, 1, 1, "1M", true);
		evidence.add(NonReferenceReadPair.create(dp[0], dp[1], ses));
		evidence.add(NonReferenceReadPair.create(dp[1], dp[0], ses));
		
		Collections.sort(calls, DirectedEvidenceOrder.ByNatural);
		Collections.sort(evidence, DirectedEvidenceOrder.ByNatural);
		ArrayList<VariantContextDirectedEvidence> result = Lists.newArrayList(new SequentialEvidenceAnnotator(pc, calls.iterator(), evidence.iterator(), 20, true, null));
		assertEquals(calls.size(), result.size());
		assertEquals(0, ((VariantContextDirectedBreakpoint)result.get(0)).getBreakpointEvidenceCountReadPair(EvidenceSubset.ALL)); // low1
		assertEquals(1, ((VariantContextDirectedBreakpoint)result.get(1)).getBreakpointEvidenceCountReadPair(EvidenceSubset.ALL)); // low2
		assertEquals(1, ((VariantContextDirectedBreakpoint)result.get(2)).getBreakpointEvidenceCountReadPair(EvidenceSubset.ALL)); // high2
		assertEquals(0, ((VariantContextDirectedBreakpoint)result.get(3)).getBreakpointEvidenceCountReadPair(EvidenceSubset.ALL)); // high1
	}
}
