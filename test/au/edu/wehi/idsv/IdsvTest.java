package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.junit.Test;

import com.google.common.base.Predicate;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.io.Files;

public class IdsvTest extends IntermediateFilesTest {
	@Test
	public void test_sv_comparision_203541() throws IOException {
		File output = new File(super.testFolder.getRoot(), "203541.vcf");
		setReference(new File("C:/dev/chr12.fa"));
		createInput(new File("testdata/203541.bam"));
		String[] args = new String[] {
				"INPUT=" + input.toString(),
				"REFERENCE=" + reference.toString(),
				"OUTPUT=" + output.toString(),
				"PER_CHR=false"
		};
		new Idsv().instanceMain(args);		
		// Should have generated two breakpoints
		List<VariantContextDirectedEvidence> ass = breaks(getVcf(getCommandlineContext(false).getFileSystemContext().getAssemblyVcf(output), null));
		ass = Lists.newArrayList(Iterables.filter(ass, new Predicate<VariantContextDirectedEvidence>() {
			@Override
			public boolean apply(VariantContextDirectedEvidence arg0) {
				return arg0.isValid() && !arg0.isFiltered();
			}
		}));
		assertEquals(2, ass.size());
		assertEquals(203476, ass.get(0).getBreakendSummary().start);
		assertEquals(FWD, ass.get(0).getBreakendSummary().direction);
		assertEquals(203540, ass.get(1).getBreakendSummary().start);
		assertEquals(BWD, ass.get(1).getBreakendSummary().direction);
		Files.copy(new File("testdata/203541.bam.idsv.realign.bam"), new File(testFolder.getRoot().getAbsolutePath() + "/input.bam.idsv.working", "input.bam.idsv.realign.bam"));
		Files.copy(new File("testdata/203541.vcf.idsv.realign.bam"), new File(testFolder.getRoot().getAbsolutePath() + "/203541.vcf.idsv.working", "203541.vcf.idsv.realign.bam"));
		
		new Idsv().instanceMain(args);
		ass = breaks(getVcf(output, null));
		assertEquals(2, ass.size());
		assertEquals(203476, ass.get(0).getBreakendSummary().start);
		assertEquals(FWD, ass.get(0).getBreakendSummary().direction);
		assertEquals(203540, ass.get(1).getBreakendSummary().start);
		assertEquals(BWD, ass.get(1).getBreakendSummary().direction);
	}
}
