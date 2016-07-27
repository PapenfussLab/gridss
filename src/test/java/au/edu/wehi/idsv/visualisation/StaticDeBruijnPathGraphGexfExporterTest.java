package au.edu.wehi.idsv.visualisation;

import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.FileFilter;
import java.io.IOException;

import org.apache.commons.io.filefilter.WildcardFileFilter;
import org.junit.Test;
import org.junit.experimental.categories.Category;

import au.edu.wehi.idsv.AssemblyAlgorithm;
import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.Hg19Tests;
import au.edu.wehi.idsv.IntermediateFilesTest;
import au.edu.wehi.idsv.ProcessStep;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMEvidenceSource;

import com.google.common.collect.ImmutableList;


public class StaticDeBruijnPathGraphGexfExporterTest extends IntermediateFilesTest {
	@Test
	@Category(Hg19Tests.class)
	public void positional_should_export_dot() throws IOException {
		File output = new File(super.testFolder.getRoot(), "chr12-244000.vcf");
		setReference(Hg19Tests.findHg19Reference("chr12.fa"));
		createInput(new File("src/test/resources/chr12-244000.bam"));
		ProcessingContext pc = getCommandlineContext(false);
		SAMEvidenceSource ses = new SAMEvidenceSource(pc, input, 0);
		ses.completeSteps(ProcessStep.ALL_STEPS);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), output);
		aes.ensureAssembled();
		File dir = new File(super.testFolder.getRoot(), "visualisation");
		File[] export = dir.listFiles((FileFilter)new WildcardFileFilter("*.dot"));
		// TODO: change graph write timing
		assertTrue(export != null && export.length > 0);
	}
	@Override
	public ProcessingContext getCommandlineContext(boolean perChr) {
		ProcessingContext pc = super.getCommandlineContext(perChr);
		pc.getAssemblyParameters().errorCorrection.maxBaseMismatchForCollapse = 1;
		pc.getAssemblyParameters().errorCorrection.collapseBubblesOnly = true;
		pc.getAssemblyParameters().method = AssemblyAlgorithm.Positional;
		pc.getAssemblyParameters().includeRemoteSplitReads = false;
		pc.getConfig().getVisualisation().assemblyGraph = true;
		pc.getConfig().getVisualisation().directory.mkdirs();
		return pc;
	}
}
