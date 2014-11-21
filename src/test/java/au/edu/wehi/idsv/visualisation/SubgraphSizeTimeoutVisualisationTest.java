package au.edu.wehi.idsv.visualisation;

import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMRecord;

import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.util.List;

import org.apache.commons.io.filefilter.WildcardFileFilter;
import org.junit.Test;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.IntermediateFilesTest;
import au.edu.wehi.idsv.ProcessStep;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMEvidenceSource;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;


public class SubgraphSizeTimeoutVisualisationTest extends IntermediateFilesTest {
	//@Test
	//@Ignore
	public void should_work_on_test_data_path_timeout_gexf() throws IOException {
		File output = new File(super.testFolder.getRoot(), "chr12-244000.vcf");
		setReference(new File("C:/dev/chr12.fa"));
		createInput(new File("src/test/resources/chr12-244000.bam"));
		SAMEvidenceSource ses = new SAMEvidenceSource(getCommandlineContext(), input, false);
		ses.completeSteps(ProcessStep.ALL_STEPS);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(getCommandlineContext(), ImmutableList.of(ses), output);
		aes.ensureAssembled();
		File dir = new File(new File(new File(super.testFolder.getRoot(), "visualisation"), "chr12"), "200000");
		File[] subgraph = dir.listFiles((FileFilter)new WildcardFileFilter("*.precollapse.gexf"));
		assertTrue(subgraph.length > 0);
	}	
	@Test
	public void should_export_on_path_timeout_gexf() throws IOException {
		ProcessingContext pc = getCommandlineContext(false);
		File output = new File(super.testFolder.getRoot(), "test_path_timeout_gexf.vcf");
		List<SAMRecord> list = Lists.newArrayList();
		for (int i = 1; i < 1000; i++) {
			SAMRecord r = new SAMRecord(pc.getBasicSamHeader());
			r.setReferenceIndex(0);
			r.setAlignmentStart(i);
			r.setReadBases(B("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
			r.setBaseQualities(B("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
			r.setCigarString("50M50S");
			r.setReadName(Integer.toString(i));
			r.setMappingQuality(40);
			list.add(r);
		}
		createInput(list.toArray(new SAMRecord[0]));
		
		SAMEvidenceSource ses = new SAMEvidenceSource(pc, input, false);
		ses.completeSteps(ProcessStep.ALL_STEPS);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), output);
		aes.ensureAssembled();
		File dir = new File(new File(new File(super.testFolder.getRoot(), "visualisation"), "polyA"), "0");
		File[] subgraph = dir.listFiles((FileFilter)new WildcardFileFilter("*.subgraph.gexf"));
		assertTrue(subgraph != null && subgraph.length > 0);
	}
	
	@Override
	public ProcessingContext getCommandlineContext(boolean perChr) {
		ProcessingContext pc = super.getCommandlineContext(perChr);
		pc.getAssemblyParameters().maxBaseMismatchForCollapse = 1;
		pc.getAssemblyParameters().collapseBubblesOnly = true;
		pc.getAssemblyParameters().debruijnGraphVisualisationDirectory = new File(super.testFolder.getRoot(), "visualisation");
		pc.getAssemblyParameters().visualiseTimeouts = true;
		pc.getAssemblyParameters().maxSubgraphFragmentWidth = 1;
		return pc;
	}
}
