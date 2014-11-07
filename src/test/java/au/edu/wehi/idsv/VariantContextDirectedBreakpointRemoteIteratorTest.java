package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import com.google.common.collect.ImmutableList;


public class VariantContextDirectedBreakpointRemoteIteratorTest extends TestHelper {
	@Test
	public void should_flip_local_remote_breakpoint_coordinates() {
		VariantContextDirectedEvidence breakend = AB().makeVariant();
		VariantContextDirectedBreakpoint breakpoint =  (VariantContextDirectedBreakpoint)AssemblyBuilder.incorporateRealignment(getContext(), breakend, Read(2, 200, 1));
		VariantContextDirectedBreakpointRemoteIterator it = new VariantContextDirectedBreakpointRemoteIterator(getContext(), AES(), ImmutableList.of(breakpoint).iterator());
		assertEquals(2, it.next().getBreakendSummary().referenceIndex);
	}
}
