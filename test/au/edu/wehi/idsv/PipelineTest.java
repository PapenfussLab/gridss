package au.edu.wehi.idsv;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.junit.Test;

import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;

import com.google.common.base.Predicate;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.io.Files;

public class PipelineTest extends CommandLineTest {
	@Test
	public void test_sv_comparision_203541() throws IOException {
		setReference(new File("C:/dev/chr12.fa"));
		createInput(new File("testdata/203541.bam"));
		extractEvidence(false);
		generateDirectedBreakpoints();
		Files.copy(new File("testdata/203541.bam.idsv.realign.bam"), new File(testFolder.getRoot().getAbsolutePath() + "/input.bam.idsv.working", "input.bam.idsv.realign.bam"));
		clusterEvidence();
		annotateBreakends();
		List<VariantContextDirectedBreakpoint> ass = Lists.newArrayList(Iterators.filter(getVcfBreaks(".idsv.vcf").iterator(), new Predicate<VariantContextDirectedBreakpoint>() {
			public boolean apply(VariantContextDirectedBreakpoint v) {
				return StringUtils.isNotBlank(v.getAssemblerProgram());
			}
		}));
		assertTrue(ass.size() > 0);
		assertEquals(203476, ass.get(0).getBreakendSummary().start);
		assertEquals(FWD, ass.get(0).getBreakendSummary().direction);
		
		assertTrue(Iterators.any(ass.iterator(), new Predicate<VariantContextDirectedBreakpoint>() {
			public boolean apply(VariantContextDirectedBreakpoint bp) {
				return bp.getStart() == 203476
						&& bp.hasAttribute("ASSLEN")
						&& bp.hasAttribute("ACONS")
						&& bp.getBreakendSummary().direction == FWD
						&& (bp.getBreakendSummary() instanceof BreakpointSummary && ((BreakpointSummary)bp.getBreakendSummary()).direction2 == BWD)
						&& bp.getAlternateAllele(0).getDisplayString().contains("203540");
			  }
		}));
	}
}
