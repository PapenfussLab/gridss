package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.junit.Ignore;
import org.junit.Test;

import com.google.common.base.Predicate;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.io.Files;

public class IdsvTest extends IntermediateFilesTest {
	@Test
	@Ignore("Moved to 3-pass system: need sc realignment bams before assemblies can be created")
	public void test_sv_comparision_203541() throws IOException {
		File output = new File(super.testFolder.getRoot(), "203541.vcf");
		setReference(new File("C:/dev/chr12.fa"));
		createInput(new File("src/test/resources/203541.bam"));
		String[] args = new String[] {
				"INPUT=" + input.toString(),
				"REFERENCE=" + reference.toString(),
				"OUTPUT=" + output.toString(),
				"TMP_DIR=" + super.testFolder.getRoot().toString(),
				"WORKING_DIR=" + super.testFolder.getRoot().toString(),
				"PER_CHR=false"
		};
		Idsv pass1 = new Idsv();
		pass1.instanceMain(args);
		Files.copy(new File("src/test/resources/203541.bam.idsv.realign.bam"), new File(testFolder.getRoot().getAbsolutePath() + "/input.bam.idsv.working", "input.bam.idsv.realign.bam"));
		Idsv pass2 = new Idsv();
		pass2.instanceMain(args);
		
		// Should have generated two breakpoints
		List<SAMRecordAssemblyEvidence> ass = Lists.newArrayList(new AssemblyEvidenceSource(pass2.getContext(), pass2.createSamEvidenceSources(), pass2.OUTPUT).iterator(false, false)); 
		ass = Lists.newArrayList(Iterables.filter(ass, new Predicate<AssemblyEvidence>() {
			@Override
			public boolean apply(AssemblyEvidence arg0) {
				
				return !arg0.isAssemblyFiltered();
			}
		}));
		assertTrue(Iterables.any(ass, new Predicate<SAMRecordAssemblyEvidence>() {
			@Override
			public boolean apply(SAMRecordAssemblyEvidence input) {
				return input.getBreakendSummary().start == 203476
						&& input.getBreakendSummary().direction == FWD;
			}
		}));
		assertTrue(Iterables.any(ass, new Predicate<SAMRecordAssemblyEvidence>() {
			@Override
			public boolean apply(SAMRecordAssemblyEvidence input) {
				return input.getBreakendSummary().start == 203540
						&& input.getBreakendSummary().direction == BWD;
			}
		}));
		Files.copy(new File("src/test/resources/203541.vcf.idsv.realign.bam"), new File(testFolder.getRoot().getAbsolutePath() + "/203541.vcf.idsv.working", "203541.vcf.idsv.realign.bam"));
		
		Idsv pass3 = new Idsv();
		pass3.instanceMain(args);
		ass = Lists.newArrayList(new AssemblyEvidenceSource(pass3.getContext(), pass3.createSamEvidenceSources(), pass3.OUTPUT).iterator(false, false));
		assertEquals(2, ass.size());
		assertEquals(203476, ass.get(0).getBreakendSummary().start);
		assertEquals(FWD, ass.get(0).getBreakendSummary().direction);
		assertEquals(203540, ass.get(1).getBreakendSummary().start);
		assertEquals(BWD, ass.get(1).getBreakendSummary().direction);
	}
	/**
	 * Either the annotator or the breakend realignment logic is moving the called breakend positions to the incorrect
	 * location.
	 */
	@Test
	@Ignore("Moved to 3-pass system: need sc realignment bams before assemblies can be created")
	public void novel_insert_2762cb5d343cd0f882d1d93a743c69a4_chr12_867621() throws IOException {
		File src = new File("src/test/resources/2762cb5d343cd0f882d1d93a743c69a4_chr12_867621");
		File target = new File(super.testFolder.getRoot(), "novelinsert");
		FileUtils.copyDirectory(src, target);
		setReference(new File("C:/dev/chr12.fa"));
		String[] args = new String[] {
				"INPUT=" + new File(target, "debug.bam").toString(),
				"REFERENCE=" + reference.toString(),
				"OUTPUT=" + new File(target, "breakend.vcf").toString(),
				"TMP_DIR=" + target.toString(),
				"WORKING_DIR=" + target.toString(),
				"CALL_ONLY_ASSEMBLIES=true",
				"PER_CHR=false"
		};
		Idsv firstPass = new Idsv();
		firstPass.instanceMain(args);
		// Should have generated two breakpoints
		List<SAMRecordAssemblyEvidence> ass = Lists.newArrayList(new AssemblyEvidenceSource(firstPass.getContext(), firstPass.createSamEvidenceSources(), firstPass.OUTPUT).iterator(false, false));
		assertEquals(2, ass.size());
	}
}
