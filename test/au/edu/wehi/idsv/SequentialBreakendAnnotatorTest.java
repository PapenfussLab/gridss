package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;

import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.SequentialBreakendAnnotator;
import au.edu.wehi.idsv.SoftClipEvidence;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.VariantContextDirectedBreakpointBuilder;
import au.edu.wehi.idsv.vcf.VcfAttributes;

import com.google.common.collect.Iterators;



public class SequentialBreakendAnnotatorTest extends TestHelper {
	public VariantContextDirectedBreakpoint go(List<DirectedEvidence> evidence, VariantContextDirectedBreakpoint toAnnotate) {
		return new SequentialBreakendAnnotator(getContext(), null, Iterators.peekingIterator(evidence.iterator())).annotate(toAnnotate);
	}
	@Test
	public void should_add_breakpoint_supporting_evidence() {
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(getContext());
		builder
			.loc("polyA", 1, 1)
			.alleles("A", "A[polyA:10[");
		VariantContextDirectedBreakpoint result = go(L(
				(DirectedEvidence)new SoftClipEvidence(getContext(), FWD, Read(0, 1, "1M3S"), Read(0, 10, "3M"))
			), builder.make());
		assertEquals(1, result.getAttributeAsInt(VcfAttributes.SOFT_CLIP_READ_COUNT.attribute(), 0));
	}
	@Test
	public void should_add_breakend_supporting_evidence() {
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(getContext());
		builder
			.loc("polyA", 1, 1)
			.alleles("A", "A[polyA:10[");
		VariantContextDirectedBreakpoint result = go(L(
				(DirectedEvidence)new SoftClipEvidence(getContext(), FWD, Read(0, 1, "1M3S"))
			), builder.make());
		assertEquals(1, result.getAttributeAsInt(VcfAttributes.SOFT_CLIP_READ_COUNT.attribute(), 0));
	}
	@Test
	public void should_use_assembly_sequence() {
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(getContext());
		builder
			.loc("polyA", 1, 1)
			.alleles("A", "A[polyA:10[");
		VariantContextDirectedBreakpoint result = go(L(
				(DirectedEvidence)AE(new BreakpointSummary(0, FWD, 1, 1, 0, BWD, 10, 10, null), 1, 2, 3, "TT")
			), builder.make());
		assertEquals("ATT[polyA:10[", result.getAlternateAllele(0).getDisplayString());
	}
	@Test
	public void should_merge_supporting_evidence() {
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(getContext());
		builder
			.loc("polyA", 1, 1)
			.alleles("A", "A[polyA:10[");
		VariantContextDirectedBreakpoint result = go(L(
				(DirectedEvidence)new SoftClipEvidence(getContext(), FWD, Read(0, 1, "1M3S"), Read(0, 10, "3M")),
				(DirectedEvidence)new SoftClipEvidence(getContext(), FWD, Read(0, 1, "1M3S"), Read(0, 10, "3M")),
				(DirectedEvidence)AE(new BreakpointSummary(0, FWD, 1, 1, 0, BWD, 10, 10, null), 1, 2, 3, "TT")
			), builder.make());
		assertEquals(2, result.getAttributeAsInt(VcfAttributes.SOFT_CLIP_READ_COUNT.attribute(), 0));
		assertEquals(1, result.getAttributeAsInt(VcfAttributes.ASSEMBLY_READS.attribute(), 0));
		assertEquals(3+3+3, result.getAttributeAsInt(VcfAttributes.REALIGN_TOTAL_LENGTH.attribute(), 0));
	}
	@Test
	public void should_not_add_breakend_nonsupporting_evidence() {
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(getContext());
		builder
			.loc("polyA", 1, 1)
			.alleles("A", "A[polyA:10[");
		VariantContextDirectedBreakpoint result = go(L(
				(DirectedEvidence)new SoftClipEvidence(getContext(), BWD, Read(0, 1, "3S1M"), Read(0, 10, "3M"))
			), builder.make());
		assertEquals(0, result.getAttributeAsInt(VcfAttributes.SOFT_CLIP_READ_COUNT.attribute(), 0));
	}
	@Test
	public void should_not_add_breakpoint_nonsupporting_evidence() {
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(getContext());
		builder
			.loc("polyA", 1, 1)
			.alleles("A", "A[polyA:10[");
		VariantContextDirectedBreakpoint result = go(L(
				(DirectedEvidence)new SoftClipEvidence(getContext(), FWD, withSequence("TTTTTTTT", Read(0, 1, "3S1M3S"))[0], withSequence("TTT", Read(0, 12, "3M"))[0])
			), builder.make());
		assertEquals(0, result.getAttributeAsInt(VcfAttributes.SOFT_CLIP_READ_COUNT.attribute(), 0));
	}
}
