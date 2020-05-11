package gridss;

import au.edu.wehi.idsv.*;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import au.edu.wehi.idsv.vcf.VcfInfoAttributes;
import com.google.common.collect.ImmutableList;
import org.junit.Test;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import static org.junit.Assert.assertEquals;

public class AnnotateInexactHomologyTest extends TestHelper {
	@Test
	public void should_calculate_inexact_homology() {
		ProcessingContext pc = getContext();
		StructuralVariationCallBuilder builder = new StructuralVariationCallBuilder(pc, new CalledBreakpointPositionLookup(), (VariantContextDirectedEvidence)new IdsvVariantContextBuilder(getContext()) {{
			breakpoint(new BreakpointSummary(2, FWD, 78, 6, BWD, 79), "");
			phredScore(50);
		}}.make());
		builder.addEvidence(SR(Read(2, 78, "1M1S"), Read(6, 79, "1M")));
		VariantContextDirectedBreakpoint e = (VariantContextDirectedBreakpoint)builder.make();
		AnnotateInexactHomology aih = new AnnotateInexactHomology();
		aih.setContext(pc);
		ExecutorService threadpool = Executors.newSingleThreadExecutor();
		e = (VariantContextDirectedBreakpoint) aih.iterator(new AutoClosingIterator<>(ImmutableList.of(e).iterator()), threadpool).next();
		assertEquals(-78, ((int[])e.getAttribute(VcfInfoAttributes.INEXACT_HOMPOS.attribute()))[0]);
		assertEquals(300, ((int[])e.getAttribute(VcfInfoAttributes.INEXACT_HOMPOS.attribute()))[1]);
		threadpool.shutdown();
	}
}
