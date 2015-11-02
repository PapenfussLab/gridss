package au.edu.wehi.idsv.debruijn;

import static org.junit.Assert.assertTrue;

import java.io.File;

import org.junit.Test;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.subgraph.DeBruijnReadGraph;
import au.edu.wehi.idsv.visualisation.DeBruijnSubgraphGexfExporter;

public class DeBruijnGraphBaseTest extends TestHelper {
	@Test
	public void should_write_gexf_() {
		AssemblyEvidenceSource aes = AES();
		aes.getContext().getConfig().getAssembly().errorCorrection.maxBaseMismatchForCollapse = 0;
		aes.getContext().getConfig().getAssembly().subgraph.traversalMaximumBranchingFactor = 10;
		aes.getContext().getConfig().getAssembly().k = 5;
		DeBruijnReadGraph g = new DeBruijnReadGraph(aes, 0, null);
		g.setGraphExporter(new DeBruijnSubgraphGexfExporter(5));
		g.addEvidence(SCE(FWD, withSequence("TAAAGTC", Read(0, 10, "4M3S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAGTCT", Read(0, 11, "3M4S"))));
		g.assembleContigsBefore(1000);
		g.addEvidence(SCE(FWD, withSequence("GAAAGTC", Read(0, 1000, "4M3S"))));
		g.addEvidence(SCE(FWD, withSequence("GAAGTCT", Read(0, 1001, "3M4S"))));
		g.addEvidence(SCE(FWD, withSequence("TAAAGTC", Read(0, 1000, "4M3S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAGTCT", Read(0, 1001, "3M4S"))));
		g.assembleContigsBefore(2000);
		File f = new File("test.gexf");
		g.getGraphExporter().saveTo(f);
		assertTrue(f.exists());
		assertTrue(f.length() > 1); // there exists SOME content in the file
		new File("test.gexf").delete();
	}
}
