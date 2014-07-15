package au.edu.wehi.idsv.debruijn.subgraph;

import static org.junit.Assert.assertEquals;

import java.util.LinkedList;
import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.AssemblyMethod;
import au.edu.wehi.idsv.AssemblyParameters;
import au.edu.wehi.idsv.AssemblyParameters.ContigAssemblyOrder;
import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;


public class SubgraphPathContigAssemblerTest extends TestHelper {

	@Test
	public void anchor_should_not_diverge_from_reference() {
		AssemblyParameters ap = new AssemblyParameters();
		ap.k = 4;
		ap.method = AssemblyMethod.DEBRUIJN_SUBGRAPH;
		ap.assemblyOrder = ContigAssemblyOrder.GreedyMaxKmer;
		ap.maxContigsPerAssembly = 1;
		
		DeBruijnReadGraph g = G(ap.k, FWD);
		g.addEvidence(SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))));
		g.addEvidence(SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))));
		g.addEvidence(NRRP(withSequence("GGTTAACC", DP(0, 1, "8M", true, 1, 10, "8M", false))));
		g.addEvidence(SCE(FWD, withSequence("TTGGT", Read(0, 10, "4M1S"))));
		//     TTAACCGGCCAATT
		//   GGTTAACC
		// TTGGT
		// ^   ^^^^ <- starts of reference kmers
		SubgraphPathContigAssembler spca = new SubgraphPathContigAssembler(g, KmerEncodingHelper.picardBaseToEncoded(ap.k, B("AATT")));
		List<LinkedList<Long>> result = spca.assembleContigs(ap);
		assertEquals("TTAACCGGCCAATT", S(g, result.get(0)));
	}
	@Test
	public void should_assemble_excluding_used_non_reference_kmers() {
		AssemblyParameters ap = new AssemblyParameters();
		ap.k = 4;
		ap.method = AssemblyMethod.DEBRUIJN_SUBGRAPH;
		ap.assemblyOrder = ContigAssemblyOrder.GreedyMaxKmer;
		ap.maxContigsPerAssembly = 1000;
		
		DeBruijnReadGraph g = G(ap.k, FWD);
		g.addEvidence(SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))));
		g.addEvidence(SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))));
		g.addEvidence(SCE(FWD, withSequence("TTAACCGAGTCCTG", Read(0, 10, "7M7S"))));
		//     TTAACCGGCCAATT
		//     TTAACCGAGTCCTG
		// TTGGT
		// ^   ^^^^ <- starts of reference kmers
		// --------^^^^^^^^^^ assembly 1 ('-' positions are trimmed)  
		// ^^^^----           assembly 2
		SubgraphPathContigAssembler spca = new SubgraphPathContigAssembler(g, KmerEncodingHelper.picardBaseToEncoded(ap.k, B("AATT")));
		List<LinkedList<Long>> result = spca.assembleContigs(ap);
		assertEquals(2, result.size());
		assertEquals("TTAACCGGCCAATT", S(g, result.get(0)));
		assertEquals("TTAACCGAGTCCTG", S(g, result.get(1)));
	}
}
