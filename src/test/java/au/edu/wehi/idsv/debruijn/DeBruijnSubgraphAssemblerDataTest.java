package au.edu.wehi.idsv.debruijn;

import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.util.CloseableIterator;

import java.io.File;
import java.util.ArrayList;

import org.junit.Test;

import au.edu.wehi.idsv.AssemblyParameters;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.IntermediateFilesTest;
import au.edu.wehi.idsv.ProcessStep;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.ReadPairParameters;
import au.edu.wehi.idsv.RealignmentParameters;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.SoftClipParameters;
import au.edu.wehi.idsv.VariantCallingParameters;
import au.edu.wehi.idsv.debruijn.subgraph.DeBruijnSubgraphAssembler;

/**
 * Tests against external data sets.
 * @author cameron.d
 *
 */
public class DeBruijnSubgraphAssemblerDataTest extends IntermediateFilesTest {
	@Test//(timeout=60000)
	public void should_assemble_MT_efficiently() {
		File hg19decoy = new File("C:\\dev\\hs37d5.fa");
		ProcessingContext pc = new ProcessingContext(
				new FileSystemContext(testFolder.getRoot(), 500000),
				new ArrayList<Header>(),
				new SoftClipParameters(),
				new ReadPairParameters(),
				new AssemblyParameters() {{
					debruijnGraphVisualisationDirectory = new File(testFolder.getRoot(), "visualisation");
					//visualiseAll = true;
					trackAlgorithmProgress = true;
					}},
				new RealignmentParameters(),
				new VariantCallingParameters(),
				hg19decoy,
				true, false);
		DeBruijnSubgraphAssembler ass = new DeBruijnSubgraphAssembler(pc, AES(pc));
		
		SAMEvidenceSource ses = new SAMEvidenceSource(new ProcessingContext(
				new FileSystemContext(new File("W:\\na12878"), 500000),
				new ArrayList<Header>(),
				new SoftClipParameters(),
				new ReadPairParameters(),
				new AssemblyParameters(),
				new RealignmentParameters(),
				new VariantCallingParameters(),
				hg19decoy,
				true, false), new File("W:\\na12878\\NA12878D_HiSeqX_R1.bam"), 0);
		ses.completeSteps(ProcessStep.ALL_STEPS); 
		CloseableIterator<DirectedEvidence> it = ses.iterator(true,  true, true, "MT");
		
		while (it.hasNext()) {
			ass.addEvidence(it.next());
		}
	}
}
