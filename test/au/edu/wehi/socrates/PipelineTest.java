package au.edu.wehi.socrates;

import static org.junit.Assert.*;
import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.junit.Test;

import com.google.common.base.Predicate;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.io.Files;
import com.sun.source.tree.AssertTree;

public class PipelineTest extends CommandLineTest {
	@Test
	public void test_sv_comparision_203541() throws IOException {
		setReference(new File("W:/refdata/genomes/chr12.fa"));
		createInput(new File("testdata/203541.bam"));
		extractEvidence(false);
		generateDirectedBreakpoints();
		Files.copy(new File("testdata/203541.bam.socrates.realign.bam"), new File(testFolder.getRoot().getAbsolutePath() + "/input.bam.socrates.working", "input.bam.socrates.realign.bam"));
		clusterEvidence();
		annotateBreakends();
		List<VariantContextDirectedBreakpoint> ass = Lists.newArrayList(Iterators.filter(getVcfBreaks(".socrates.vcf").iterator(), new Predicate<VariantContextDirectedBreakpoint>() {
			public boolean apply(VariantContextDirectedBreakpoint v) {
				return StringUtils.isNotBlank(v.getAssemblerProgram());
			}
		}));
		assertTrue(ass.size() > 0);
		assertEquals(203476, ass.get(0).getBreakendSummary().start);
		assertEquals(FWD, ass.get(0).getBreakendSummary().direction);
	}
}
