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
	public void should_add_breakend_supporting_evidence() {
		VariantContextDirectedBreakpointBuilder builder = new VariantContextDirectedBreakpointBuilder(getContext());
		builder
			.chr("polyA")
			.start(1)
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
			.chr("polyA")
			.start(1)
			.alleles("A", "A[polyA:10[");
		VariantContextDirectedBreakpoint result = go(L(
				(DirectedEvidence)AE(new BreakpointSummary(0, FWD, 1, 1, 0, BWD, 10, 10, null), 1, 2, 3, "TT")
			), builder.make());
		assertEquals("ATT[polyA:10[", result.getAlternateAllele(0).getDisplayString());
	}
}
