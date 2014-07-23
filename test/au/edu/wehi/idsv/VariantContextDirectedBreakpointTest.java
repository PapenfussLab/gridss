package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import au.edu.wehi.idsv.vcf.VcfAttributes;


public class VariantContextDirectedBreakpointTest extends TestHelper {
	@Test
	public void getAnchorSequenceString_should_return_entire_assembly_anchor() {
		assertEquals("GAA", ((VariantContextDirectedBreakpoint)minimalBreakend()
				.breakpoint(new BreakpointSummary(0, FWD, 1, 1, 1, FWD, 5, 5), "T")
				.attribute(VcfAttributes.ASSEMBLY_CONSENSUS.attribute(), "GAATT")
				.attribute(VcfAttributes.ASSEMBLY_LENGTH_LOCAL_MAX.attribute(), 3)
				.make())
			.getAnchorSequenceString());
		assertEquals("ATT", ((VariantContextDirectedBreakpoint)minimalBreakend()
				.breakpoint(new BreakpointSummary(0, BWD, 1, 1, 1, FWD, 5, 5), "")
				.attribute(VcfAttributes.ASSEMBLY_CONSENSUS.attribute(), "GAATT")
				.attribute(VcfAttributes.ASSEMBLY_LENGTH_LOCAL_MAX.attribute(), 3)
				.make())
			.getAnchorSequenceString());
	}
}
