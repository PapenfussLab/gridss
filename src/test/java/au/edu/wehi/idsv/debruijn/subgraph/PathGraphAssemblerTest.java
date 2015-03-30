package au.edu.wehi.idsv.debruijn.subgraph;

import static org.junit.Assert.assertEquals;

import java.util.LinkedList;
import java.util.List;

import org.junit.Ignore;
import org.junit.Test;

import au.edu.wehi.idsv.AssemblyMethod;
import au.edu.wehi.idsv.AssemblyParameters;
import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.visualisation.NontrackingSubgraphTracker;


public class PathGraphAssemblerTest extends TestHelper {
	@Test
	public void anchor_should_not_diverge_from_reference() {
		AssemblyParameters ap = new AssemblyParameters();
		ap.k = 4;
		ap.method = AssemblyMethod.DEBRUIJN_SUBGRAPH;
		ap.subgraphAssemblyTraversalMaximumBranchingFactor = 1;
		ap.maxContigsPerAssembly = 1;
		ap.assemblyBackToReference = false;
		
		DeBruijnReadGraph g = G(ap.k, FWD);
		g.addEvidence(SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))));
		g.addEvidence(SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))));
		g.addEvidence(NRRP(withSequence("GGTTAACC", DP(0, 1, "8M", true, 1, 10, "8M", false))));
		g.addEvidence(SCE(FWD, withSequence("TTGGT", Read(0, 10, "4M1S"))));
		//     TTAACCGGCCAATT
		//   GGTTAACC
		// TTGGT
		// ^   ^^^^ <- starts of reference kmers
		PathGraphAssembler pga = new PathGraphAssembler(g, ap, KmerEncodingHelper.picardBaseToEncoded(ap.k, B("AATT")), new NontrackingSubgraphTracker());
		List<LinkedList<Long>> result = pga.assembleContigs();
		assertEquals(2, result.size());
		assertEquals("TTAACCGGCCAATT", S(g, result.get(0)));
	}
	@Ignore("now removing entire non-reference subgraph after each assembly")
	//@Test
	public void should_assemble_excluding_used_non_reference_kmers() {
		AssemblyParameters ap = new AssemblyParameters();
		ap.k = 4;
		ap.method = AssemblyMethod.DEBRUIJN_SUBGRAPH;
		ap.subgraphAssemblyTraversalMaximumBranchingFactor = 1;
		ap.maxContigsPerAssembly = 1000;
		ap.maxBaseMismatchForCollapse = 0;
		ap.assemblyBackToReference = false;
		
		DeBruijnReadGraph g = G(ap.k, FWD);
		g.addEvidence(SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))));
		g.addEvidence(SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))));
		g.addEvidence(NRRP(withSequence(           "GCCAATC", OEA(0, 10, "7M", true)))); // not on main patch, but reachable
		g.addEvidence(SCE(FWD, withSequence("TTAACCGAGTCCTG", Read(0, 10, "7M7S"))));
		
		PathGraphAssembler pga = new PathGraphAssembler(g, ap, KmerEncodingHelper.picardBaseToEncoded(ap.k, B("AATT")), new NontrackingSubgraphTracker());
		List<LinkedList<Long>> result = pga.assembleContigs();
		assertEquals(3, result.size());
		assertEquals("TTAACCGGCCAATT", S(g, result.get(0)));
		assertEquals("TTAACCGAGTCCTG", S(g, result.get(1)));
	}
	@Test
	public void should_assemble_excluding_all_reachable_non_reference_kmers() {
		AssemblyParameters ap = new AssemblyParameters();
		ap.k = 4;
		ap.method = AssemblyMethod.DEBRUIJN_SUBGRAPH;
		ap.subgraphAssemblyTraversalMaximumBranchingFactor = 1;
		ap.maxContigsPerAssembly = 1000;
		ap.maxBaseMismatchForCollapse = 0;
		ap.assemblyBackToReference = false;
		
		DeBruijnReadGraph g = G(ap.k, FWD);
		g.addEvidence(SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))));
		g.addEvidence(SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))));
		g.addEvidence(NRRP(withSequence(              "AATC", OEA(0, 10, "7M", true)))); // not on main patch, but reachable
		g.addEvidence(SCE(FWD, withSequence("TTAACCGAGTCCTG", Read(0, 10, "7M7S"))));
		
		PathGraphAssembler pga = new PathGraphAssembler(g, ap, KmerEncodingHelper.picardBaseToEncoded(ap.k, B("AATT")), new NontrackingSubgraphTracker());
		List<LinkedList<Long>> result = pga.assembleContigs();
		assertEquals(2, result.size());
		assertEquals("TTAACCGGCCAATT", S(g, result.get(0)));
		assertEquals("TTAACCGAGTCCTG", S(g, result.get(1)));
	}
	@Test
	public void should_maximise_non_reference_weight() {
		AssemblyParameters ap = new AssemblyParameters();
		ap.k = 4;
		ap.method = AssemblyMethod.DEBRUIJN_SUBGRAPH;
		ap.maxContigsPerAssembly = 1000;
		ap.maxBaseMismatchForCollapse = 0;
		ap.assemblyBackToReference = false;
		DeBruijnReadGraph g = G(ap.k, FWD);
		g.addEvidence(SCE(FWD, withSequence("AAATGGGG", Read(0, 10, "6M3S"))));
		g.addEvidence(SCE(FWD, withSequence("GTACCCGGGG", Read(0, 10, "9M2S"))));
		// should take the shorter assembly that has a longer soft clip
		PathGraphAssembler pga = new PathGraphAssembler(g, ap, KmerEncodingHelper.picardBaseToEncoded(ap.k, B("GGGG")), new NontrackingSubgraphTracker());
		List<LinkedList<Long>> result = pga.assembleContigs();
		assertEquals(1, result.size());
		assertEquals("AAATGGGG", S(g, result.get(0)));
	}
	@Test
	public void should_maximise_reference_weight() {
		AssemblyParameters ap = new AssemblyParameters();
		ap.k = 4;
		ap.method = AssemblyMethod.DEBRUIJN_SUBGRAPH;
		ap.subgraphAssemblyTraversalMaximumBranchingFactor = 1;
		ap.maxContigsPerAssembly = 1000;
		ap.maxBaseMismatchForCollapse = 0;
		ap.assemblyBackToReference = false;
		DeBruijnReadGraph g = G(ap.k, FWD);
		g.addEvidence(SCE(FWD, withSequence( "TAAATGGG", Read(0, 10, "7M1S"))));
		g.addEvidence(SCE(FWD, withSequence( "TAAATGGG", Read(0, 10, "7M1S"))));
		g.addEvidence(SCE(FWD, withSequence( "TAAATGGG", Read(0, 10, "7M1S"))));
		g.addEvidence(SCE(FWD, withSequence("GCAAATGGG", Read(0, 10, "8M1S"))));
		g.addEvidence(SCE(FWD, withSequence("GCAAATGGG", Read(0, 10, "8M1S"))));
		g.addEvidence(SCE(FWD, withSequence("GCAAATGGG", Read(0, 10, "8M1S"))));
		g.addEvidence(SCE(FWD, withSequence("CCAAATGGG", Read(0, 10, "8M1S")))); // worst weight at the immediate next kmer
		g.addEvidence(SCE(FWD, withSequence("GGTACCCAAATGGG", Read(0, 10, "13M1S")))); // but best overall
		// should take the shorter assembly that has a longer soft clip
		PathGraphAssembler pga = new PathGraphAssembler(g, ap, KmerEncodingHelper.picardBaseToEncoded(ap.k, B("TGGG")), new NontrackingSubgraphTracker());
		List<LinkedList<Long>> result = pga.assembleContigs();
		assertEquals(1, result.size());
		assertEquals("GGTACCCAAATGGG", S(g, result.get(0)));
	}
	@Test
	public void should_assemble_simple_sequence() {
		AssemblyParameters ap = new AssemblyParameters();
		ap.k = 4;
		ap.method = AssemblyMethod.DEBRUIJN_SUBGRAPH;
		ap.subgraphAssemblyTraversalMaximumBranchingFactor = 16;
		ap.maxContigsPerAssembly = 1000;
		ap.maxBaseMismatchForCollapse = 0;
		ap.assemblyBackToReference = false;
		DeBruijnReadGraph g = G(ap.k, FWD);
		g.addEvidence(SCE(FWD, withSequence("GGTACCCAAATGGG", Read(0, 10, "10M4S"))));
		PathGraphAssembler pga = new PathGraphAssembler(g, ap, KmerEncodingHelper.picardBaseToEncoded(ap.k, B("TGGG")), new NontrackingSubgraphTracker());
		List<LinkedList<Long>> result = pga.assembleContigs();
		assertEquals(1, result.size());
		assertEquals("GGTACCCAAATGGG", S(g, result.get(0)));
	}
	/**
	 * Assembly back into reference is required for small-medium insertion detection
	 */
	@Test
	public void assemblyBackToReference_should_control_assembly_from_breakend_back_to_reference_kmers() {
		AssemblyParameters ap = new AssemblyParameters();
		ap.k = 3;
		ap.method = AssemblyMethod.DEBRUIJN_SUBGRAPH;
		ap.subgraphAssemblyTraversalMaximumBranchingFactor = 16;
		ap.maxContigsPerAssembly = 1000;
		ap.maxBaseMismatchForCollapse = 0;
		ap.assemblyBackToReference = false;
		DeBruijnReadGraph g = G(ap.k, FWD);
		g.addEvidence(SCE(FWD, withSequence("ACTGT", Read(0, 10, "3M2S"))));
		g.addEvidence(SCE(FWD, withSequence("GTCA", Read(0, 10, "3M1S"))));
		PathGraphAssembler pga = new PathGraphAssembler(g, ap, KmerEncodingHelper.picardBaseToEncoded(ap.k, B("ACT")), new NontrackingSubgraphTracker());
		List<LinkedList<Long>> result = pga.assembleContigs();
		assertEquals(2, result.size());
		assertEquals("ACTGT", S(g, result.get(0)));
		
		ap.assemblyBackToReference = true;
		pga = new PathGraphAssembler(g, ap, KmerEncodingHelper.picardBaseToEncoded(ap.k, B("ACT")), new NontrackingSubgraphTracker());
		result = pga.assembleContigs();
		assertEquals(2, result.size());
		assertEquals("ACTGTC", S(g, result.get(0)));
		// Anchor     ACT>
		// Breakend    CTG>
		//              TGT>
		// Ref           GTC>
		//                TCA is not included as that's back out to another breakend
	}
	@Test
	public void should_order_assembly_by_breakend_weight() {
		AssemblyParameters ap = new AssemblyParameters();
		ap.k = 4;
		ap.method = AssemblyMethod.DEBRUIJN_SUBGRAPH;
		ap.subgraphAssemblyTraversalMaximumBranchingFactor = 1;
		ap.maxContigsPerAssembly = 4;
		ap.maxBaseMismatchForCollapse = 0;
		ap.assemblyBackToReference = false;
		DeBruijnReadGraph g = G(ap.k, FWD);
		//                                    RRRRRR
		g.addEvidence(SCE(FWD, withSequence( "TTTTCGAGAT", Read(0, 10, "8M1S"))));
		g.addEvidence(SCE(FWD, withSequence( "TTTTCGCCGGT", Read(0, 10, "8M3S"))));
		g.addEvidence(SCE(FWD, withSequence( "TTTTCGTTAACC", Read(0, 10, "8M6S"))));
		PathGraphAssembler pga = new PathGraphAssembler(g, ap, KmerEncodingHelper.picardBaseToEncoded(ap.k, B("TTTT")), new au.edu.wehi.idsv.visualisation.NontrackingSubgraphTracker());
		List<LinkedList<Long>> result = pga.assembleContigs();
		assertEquals(3, result.size());
		assertEquals("TTTTCGTTAACC", S(g, result.get(0)));
		assertEquals("TTTTCGCCGGT", S(g, result.get(1)));
		assertEquals("TTTTCGAGAT", S(g, result.get(2)));
	}
	@Test
	public void should_fall_back_to_greedy_assemble() {
		AssemblyParameters ap = new AssemblyParameters();
		ap.k = 4;
		ap.method = AssemblyMethod.DEBRUIJN_SUBGRAPH;
		ap.maxContigsPerAssembly = 1000;
		ap.maxBaseMismatchForCollapse = 0;
		ap.maxPathTraversalNodes = 2;
		ap.assemblyBackToReference = false;
		DeBruijnReadGraph g = G(ap.k, FWD);
		g.addEvidence(SCE(FWD, withSequence("AAATGGGGGGGGGGG", Read(0, 10, "6M10S"))));
		g.addEvidence(SCE(FWD, withSequence("AAATGGGC", Read(0, 10, "6M3S"))));
		g.addEvidence(SCE(FWD, withSequence("AAATGGGTAA", Read(0, 10, "6M5S"))));
		g.addEvidence(SCE(FWD, withSequence("AAATGGGTCC", Read(0, 10, "6M5S"))));
		g.addEvidence(NRRP(withSequence("GGAC", OEA(0, 10, "4M", true)))); // GTCC reverse comp GGAC
		g.addEvidence(NRRP(withSequence("GGAC", OEA(0, 10, "4M", true))));
		g.addEvidence(SCE(FWD, withSequence("AAATGGGA", Read(0, 10, "6M3S"))));
		g.addEvidence(SCE(FWD, withSequence("AAATGGGA", Read(0, 10, "6M3S"))));
		g.addEvidence(SCE(FWD, withSequence("AAATGGGA", Read(0, 10, "6M3S"))));
		g.addEvidence(SCE(FWD, withSequence("AAATGGGA", Read(0, 10, "6M3S"))));
		// greedy traversal
		PathGraphAssembler pga = new PathGraphAssembler(g, ap, KmerEncodingHelper.picardBaseToEncoded(ap.k, B("GGGG")), new NontrackingSubgraphTracker());
		//new StaticDeBruijnSubgraphPathGraphGexfExporter(ap.k)
		//	.snapshot(pga)
		//	.saveTo(new File("C:\\temp\\test.gexf"));
		List<LinkedList<Long>> result = pga.assembleContigs();
		assertEquals(1, result.size());
		assertEquals("AAATGGGGA", S(g, result.get(0)));
	}
}