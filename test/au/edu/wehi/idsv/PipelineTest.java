package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.junit.Test;

import com.google.common.base.Predicate;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.io.Files;

public class PipelineTest extends IntermediateFilesTest {
	@Test
	public void test_sv_comparision_203541() throws IOException {
		setReference(new File("C:/dev/chr12.fa"));
		createInput(new File("testdata/203541.bam"));
		SAMEvidenceSource ses = new SAMEvidenceSource(getContext(), input, false);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(getCommandlineContext(), ImmutableList.of(ses), new File(super.testFolder.getRoot(), "out.vcf"));
		aes.ensureAssembled();
		// Should have generated two breakpoints
		List<VariantContextDirectedBreakpoint> ass = breaks(getAssembly(aes));
		assertTrue(ass.size() == 2);
		assertEquals(203476, ass.get(0).getBreakendSummary().start);
		assertEquals(FWD, ass.get(0).getBreakendSummary().direction);
		assertEquals(203540, ass.get(1).getBreakendSummary().start);
		assertEquals(BWD, ass.get(1).getBreakendSummary().direction);
		Files.copy(new File("testdata/203541.bam.idsv.realign.bam"), new File(testFolder.getRoot().getAbsolutePath() + "/input.bam.idsv.working", "input.bam.idsv.realign.bam"));
		
		// TODO: follow up pipeline
		clusterEvidence();
		annotateBreakends();
		ass = Lists.newArrayList(Iterators.filter(getVcfBreaks(".idsv.vcf").iterator(), new Predicate<VariantContextDirectedBreakpoint>() {
			public boolean apply(VariantContextDirectedBreakpoint v) {
				return StringUtils.isNotBlank(v.getAssemblerProgram());
			}
		}));
		assertTrue(ass.size() == 2);
		assertEquals(203476, ass.get(0).getBreakendSummary().start);
		assertEquals(FWD, ass.get(0).getBreakendSummary().direction);
		assertEquals(203540, ass.get(1).getBreakendSummary().start);
		assertEquals(BWD, ass.get(1).getBreakendSummary().direction);
	}
}
