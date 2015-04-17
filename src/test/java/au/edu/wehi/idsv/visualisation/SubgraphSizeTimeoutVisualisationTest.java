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
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;


public class SubgraphSizeTimeoutVisualisationTest extends IntermediateFilesTest {
	//@Test
	//@Ignore
	public void should_work_on_test_data_path_timeout_gexf() throws IOException {
		File output = new File(super.testFolder.getRoot(), "chr12-244000.vcf");
		setReference(new File("C:/dev/chr12.fa"));
		createInput(new File("src/test/resources/chr12-244000.bam"));
		SAMEvidenceSource ses = new SAMEvidenceSource(getCommandlineContext(), input, 0);
		ses.completeSteps(ProcessStep.ALL_STEPS);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(getCommandlineContext(), ImmutableList.of(ses), output);
		aes.ensureAssembled();
		File dir = new File(new File(new File(super.testFolder.getRoot(), "visualisation"), "chr12"), "200000");
		File[] subgraph = dir.listFiles((FileFilter)new WildcardFileFilter("*.precollapse.gexf"));
		assertTrue(subgraph.length > 0);
	}	
	@Test
	public void should_export_on_path_timeout_gexf() throws IOException {
		int k = 6;
		ProcessingContext pc = getCommandlineContext(false);
		pc.getAssemblyParameters().maxPathTraversalNodes = 1024;
		pc.getAssemblyParameters().k = k;
		pc.getSoftClipParameters().minLength = 1;
		File output = new File(super.testFolder.getRoot(), "test_path_timeout_gexf.vcf");
		List<SAMRecord> list = Lists.newArrayList();
		for (int i = 0; i < 1 << (2*k); i++) {
			SAMRecord r = new SAMRecord(pc.getBasicSamHeader());
			r.setReferenceIndex(0);
			r.setAlignmentStart(1);
			r.setReadBases(KmerEncodingHelper.encodedToPicardBases(k, i));
			r.setBaseQualities(B(40, k));
			r.setCigarString(String.format("1M%dS", k-1));
			r.setReadName(Integer.toString(i));
			r.setMappingQuality(40);
			list.add(r);
		}
		createInput(list.toArray(new SAMRecord[0]));
		
		SAMEvidenceSource ses = new SAMEvidenceSource(pc, input, 0);
		ses.completeSteps(ProcessStep.ALL_STEPS);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), output);
		aes.ensureAssembled();
		assertTrue(new File(new File(super.testFolder.getRoot(), "visualisation"), "pathTraversalTimeout-1.subgraph.gexf").exists());
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
