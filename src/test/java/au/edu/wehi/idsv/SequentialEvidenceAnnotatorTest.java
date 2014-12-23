package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
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
			.loc("polyA", 1, 1)
			.alleles("A", "A[polyA:10[");
		VariantContextDirectedBreakpoint result = (VariantContextDirectedBreakpoint)go(L(
				(DirectedEvidence)SoftClipEvidence.create(getContext(), SES(), FWD, Read(0, 1, "1M3S"), Read(0, 10, "3M"))
			), (VariantContextDirectedEvidence)builder.make());
		assertEquals(1, result.getBreakpointEvidenceCountSoftClip(null));
	}
	@Test
	public void should_add_breakend_supporting_evidence() {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext());
		builder
			.loc("polyA", 1, 1)
			.alleles("A", "A[polyA:10[");
		VariantContextDirectedEvidence result = go(L(
				(DirectedEvidence)SoftClipEvidence.create(getContext(), SES(), FWD, Read(0, 1, "1M3S"))
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
			.loc("polyA", 1, 1)
			.alleles("A", "A[polyA:10[");
		VariantContextDirectedBreakpoint result = (VariantContextDirectedBreakpoint)go(L(
				(DirectedEvidence)SoftClipEvidence.create(getContext(), SES(), FWD, Read(0, 1, "1M2S"), Read(0, 10, "3M")),
				(DirectedEvidence)SoftClipEvidence.create(getContext(), SES(), FWD, Read(0, 1, "1M3S"), Read(0, 10, "3M"))
			), (VariantContextDirectedEvidence)builder.make());
		assertEquals(2, result.getBreakpointEvidenceCountSoftClip(null));
	}
	@Test
	public void should_incorporate_sc_margin_when_adding_evidence() {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext());
		builder
			.loc("polyA", 2, 2)
			.alleles("A", "A[polyA:10[");
		VariantContextDirectedBreakpoint result = (VariantContextDirectedBreakpoint)go(L(
				(DirectedEvidence)SoftClipEvidence.create(getContext(), SES(), FWD, Read(0, 1, "1M3S"))
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
		evidence.add(SoftClipEvidence.create(getContext(), SES(), FWD, withSequence("NNNNN", Read(0, 12, "1M4S"))[0], withSequence("NNN", Read(0, 10, "3M"))[0]));
		evidence.add(SoftClipEvidence.create(getContext(), SES(), FWD, withSequence("NNNN", Read(0, 12, "1M3S"))[0], withSequence("NNN", Read(0, 10, "3M"))[0]));
		
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
			.loc("polyA", 1, 1)
			.alleles("A", "A[polyA:10[");
		VariantContextDirectedBreakpoint result = (VariantContextDirectedBreakpoint)go(L(
				(DirectedEvidence)SoftClipEvidence.create(getContext(), SES(), BWD, Read(0, 1, "3S1M"), Read(0, 10, "3M"))
			), (VariantContextDirectedEvidence)builder.make());
		assertEquals(0, result.getBreakpointEvidenceCountSoftClip(null));
		assertEquals(0, result.getBreakendEvidenceCountSoftClip(null));
	}
	@Test
	public void should_not_add_breakpoint_nonsupporting_evidence() {
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(getContext());
		builder
			.loc("polyA", 1, 1)
			.alleles("A", "A[polyA:10[");
		VariantContextDirectedBreakpoint result = (VariantContextDirectedBreakpoint)go(L(
				(DirectedEvidence)SoftClipEvidence.create(getContext(), SES(), FWD, withSequence("TTTTTTTT", Read(0, 5, "3S1M3S"))[0], withSequence("TTT", Read(0, 12, "3M"))[0])
			), (VariantContextDirectedEvidence)builder.make());
		assertEquals(0, result.getBreakpointEvidenceCountSoftClip(null));
		assertEquals(0, result.getBreakendEvidenceCountSoftClip(null));
	}
}
