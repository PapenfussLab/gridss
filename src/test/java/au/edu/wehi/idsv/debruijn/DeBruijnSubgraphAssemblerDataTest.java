package au.edu.wehi.idsv.debruijn;

import java.io.File;
import java.util.ArrayList;

import org.apache.commons.configuration.ConfigurationException;
import org.junit.Ignore;
import org.junit.Test;
import org.junit.experimental.categories.Category;

import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.Hg19Tests;
import au.edu.wehi.idsv.IntermediateFilesTest;
import au.edu.wehi.idsv.ProcessStep;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.configuration.GridssConfiguration;
import au.edu.wehi.idsv.debruijn.subgraph.DeBruijnSubgraphAssembler;
import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.util.CloseableIterator;

/**
 * Tests against external data sets.
 * @author Daniel Cameron
 *
 */
public class DeBruijnSubgraphAssemblerDataTest extends IntermediateFilesTest {
	@Test//(timeout=60000)
	@Ignore() // TODO: fix testing architecture so performance test requiring external data are configurable
	@Category(Hg19Tests.class)
	public void should_assemble_MT_efficiently() throws ConfigurationException {
		File hg19decoy = Hg19Tests.findHg19Reference("hs37d5.fa");
		ProcessingContext pc = new ProcessingContext(
				new FileSystemContext(testFolder.getRoot(), 500000), hg19decoy, true, null,
				new ArrayList<Header>(), new GridssConfiguration((File)null, testFolder.getRoot()));
		pc.getConfig().getVisualisation().assemblyProgress = true;
		
		DeBruijnSubgraphAssembler ass = new DeBruijnSubgraphAssembler(AES(pc));
		
		SAMEvidenceSource ses = new SAMEvidenceSource(new ProcessingContext(
				new FileSystemContext(new File("W:\\na12878"), 500000), hg19decoy, true, null,
				new ArrayList<Header>(), new GridssConfiguration()),
			new File("W:\\na12878\\NA12878D_HiSeqX_R1.bam"), 0);
		ses.completeSteps(ProcessStep.ALL_STEPS); 
		CloseableIterator<DirectedEvidence> it = ses.iterator(true, true, true, "MT");
		
		while (it.hasNext()) {
			ass.addEvidence(it.next());
		}
	}
}
