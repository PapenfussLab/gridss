package au.edu.wehi.idsv;

import au.edu.wehi.idsv.picard.InMemoryReferenceSequenceFile;
import htsjdk.samtools.metrics.Header;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;

public class VcfBreakendSummaryTest extends TestHelper {
	@Test
	public void should_allow_colon_in_breakpoint_contig_name() {
		InMemoryReferenceSequenceFile ref = new InMemoryReferenceSequenceFile(new String[] { "HLA-A*01:04N:1-HLA-A*02:95:3388" }, new byte[][] { B("ACGT") });
		ProcessingContext pc = new ProcessingContext(getFSContext(), null, ref, new ArrayList<Header>(), getConfig());
		IdsvVariantContextBuilder builder = new IdsvVariantContextBuilder(pc);
		builder.chr("HLA-A*01:04N:1-HLA-A*02:95:3388").start(1).stop(1).alleles("A", "A.");
		builder.breakpoint(new BreakpointSummary(0, FWD, 1, 1, 1, 0, BWD, 2, 2, 2), "");
		IdsvVariantContext vc = builder.make();
		VcfBreakendSummary summary = new VcfBreakendSummary(pc, vc);
		Assert.assertEquals(0, summary.location.referenceIndex);
	}
}
