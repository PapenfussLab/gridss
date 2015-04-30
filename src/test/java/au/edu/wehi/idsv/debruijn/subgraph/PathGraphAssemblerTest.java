package au.edu.wehi.idsv.debruijn.subgraph;

import static org.junit.Assert.assertEquals;

import java.util.List;

import org.junit.Ignore;
import org.junit.Test;

import au.edu.wehi.idsv.AssemblyParameters;
import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.debruijn.DeBruijnPathNode;
import au.edu.wehi.idsv.debruijn.DeBruijnPathNodeFactory;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.graph.PathNode;
import au.edu.wehi.idsv.visualisation.NontrackingSubgraphTracker;

import com.google.common.base.Function;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;


public class PathGraphAssemblerTest extends TestHelper {
	PathGraphAssembler<DeBruijnSubgraphNode, DeBruijnPathNode<DeBruijnSubgraphNode>> PGA(DeBruijnReadGraph g, AssemblyParameters ap, String seed) {
		DeBruijnSubgraphNode seedNode = g.getKmer(KmerEncodingHelper.picardBaseToEncoded(ap.k, B(seed)));
		PathGraphAssembler<DeBruijnSubgraphNode, DeBruijnPathNode<DeBruijnSubgraphNode>> pga = new PathGraphAssembler<DeBruijnSubgraphNode, DeBruijnPathNode<DeBruijnSubgraphNode>>(
			g,
			ap,
			ImmutableList.of(seedNode),
			new DeBruijnPathNodeFactory<DeBruijnSubgraphNode>(g),
			new NontrackingSubgraphTracker<DeBruijnSubgraphNode, DeBruijnPathNode<DeBruijnSubgraphNode>>());
		return pga;
	}
	private List<List<DeBruijnSubgraphNode>> flatten(List<List<DeBruijnPathNode<DeBruijnSubgraphNode>>> list) {
		List<List<DeBruijnSubgraphNode>> result = Lists.newArrayList(Iterables.transform(list, new Function<List<DeBruijnPathNode<DeBruijnSubgraphNode>>, List<DeBruijnSubgraphNode>>() {
			@Override
			public List<DeBruijnSubgraphNode> apply(List<DeBruijnPathNode<DeBruijnSubgraphNode>> input) {
				return Lists.newArrayList(PathNode.nodeIterator(input.iterator()));
			}
		}));
		return result;
	}
	@Test
	public void anchor_should_not_diverge_from_reference() {
		AssemblyParameters ap = new AssemblyParameters();
		ap.k = 4;
		ap.subgraphAssemblyTraversalMaximumBranchingFactor = 1;
		ap.maxContigsPerAssembly = 1;
		
		DeBruijnReadGraph g = RG(ap.k);
		g.addEvidence(SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))));
		g.addEvidence(SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))));
		g.addEvidence(NRRP(withSequence("GGTTAACC", DP(0, 1, "8M", true, 1, 10, "8M", false))));
		g.addEvidence(SCE(FWD, withSequence("TTGGT", Read(0, 10, "4M1S"))));
		//     TTAACCGGCCAATT
		//   GGTTAACC
		// TTGGT
		// ^   ^^^^ <- starts of reference kmers
		PathGraphAssembler<DeBruijnSubgraphNode, DeBruijnPathNode<DeBruijnSubgraphNode>> pga = PGA(g, ap, "AATT");
		List<List<DeBruijnSubgraphNode>> result = flatten(pga.assembleContigs());
		assertEquals(2, result.size());
		assertEquals("TTAACCGGCCAATT", S(g, result.get(0)));
	}
	@Ignore("now removing entire non-reference subgraph after each assembly")
	//@Test
	public void should_assemble_excluding_used_non_reference_kmers() {
		AssemblyParameters ap = new AssemblyParameters();
		ap.k = 4;
		ap.subgraphAssemblyTraversalMaximumBranchingFactor = 1;
		ap.maxContigsPerAssembly = 1000;
		ap.maxBaseMismatchForCollapse = 0;
		
		DeBruijnReadGraph g = RG(ap.k);
		g.addEvidence(SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))));
		g.addEvidence(SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))));
		g.addEvidence(NRRP(withSequence(           "GCCAATC", OEA(0, 10, "7M", true)))); // not on main patch, but reachable
		g.addEvidence(SCE(FWD, withSequence("TTAACCGAGTCCTG", Read(0, 10, "7M7S"))));
		
		PathGraphAssembler<DeBruijnSubgraphNode, DeBruijnPathNode<DeBruijnSubgraphNode>> pga = PGA(g, ap, "AATT");
		List<List<DeBruijnSubgraphNode>> result = flatten(pga.assembleContigs());
		assertEquals(3, result.size());
		assertEquals("TTAACCGGCCAATT", S(g, result.get(0)));
		assertEquals("TTAACCGAGTCCTG", S(g, result.get(1)));
	}
	@Test
	public void should_assemble_excluding_all_reachable_non_reference_kmers() {
		AssemblyParameters ap = new AssemblyParameters();
		ap.k = 4;
		ap.subgraphAssemblyTraversalMaximumBranchingFactor = 1;
		ap.maxContigsPerAssembly = 1000;
		ap.maxBaseMismatchForCollapse = 0;
		
		
		DeBruijnReadGraph g = RG(ap.k);
		g.addEvidence(SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))));
		g.addEvidence(SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))));
		g.addEvidence(NRRP(withSequence(              "AATC", OEA(0, 10, "7M", true)))); // not on main patch, but reachable
		g.addEvidence(SCE(FWD, withSequence("TTAACCGAGTCCTG", Read(0, 10, "7M7S"))));
		
		PathGraphAssembler<DeBruijnSubgraphNode, DeBruijnPathNode<DeBruijnSubgraphNode>> pga = PGA(g, ap, "AATT");
		List<List<DeBruijnSubgraphNode>> result = flatten(pga.assembleContigs());
		assertEquals(2, result.size());
		assertEquals("TTAACCGGCCAATT", S(g, result.get(0)));
		assertEquals("TTAACCGAGTCCTG", S(g, result.get(1)));
	}
	@Test
	public void should_maximise_non_reference_weight() {
		AssemblyParameters ap = new AssemblyParameters();
		ap.k = 4;
		ap.maxContigsPerAssembly = 1000;
		ap.maxBaseMismatchForCollapse = 0;
		
		DeBruijnReadGraph g = RG(ap.k);
		g.addEvidence(SCE(FWD, withSequence("AAATGGGG", Read(0, 10, "6M3S"))));
		g.addEvidence(SCE(FWD, withSequence("GTACCCGGGG", Read(0, 10, "9M2S"))));
		// should take the shorter assembly that has a longer soft clip
		PathGraphAssembler<DeBruijnSubgraphNode, DeBruijnPathNode<DeBruijnSubgraphNode>> pga = PGA(g, ap, "GGGG");
		List<List<DeBruijnSubgraphNode>> result = flatten(pga.assembleContigs());
		assertEquals(1, result.size());
		assertEquals("AAATGGGG", S(g, result.get(0)));
	}
	@Test
	public void should_maximise_reference_weight() {
		AssemblyParameters ap = new AssemblyParameters();
		ap.k = 4;
		ap.subgraphAssemblyTraversalMaximumBranchingFactor = 1;
		ap.maxContigsPerAssembly = 1000;
		ap.maxBaseMismatchForCollapse = 0;
		
		DeBruijnReadGraph g = RG(ap.k);
		g.addEvidence(SCE(FWD, withSequence( "TAAATGGG", Read(0, 10, "7M1S"))));
		g.addEvidence(SCE(FWD, withSequence( "TAAATGGG", Read(0, 10, "7M1S"))));
		g.addEvidence(SCE(FWD, withSequence( "TAAATGGG", Read(0, 10, "7M1S"))));
		g.addEvidence(SCE(FWD, withSequence("GCAAATGGG", Read(0, 10, "8M1S"))));
		g.addEvidence(SCE(FWD, withSequence("GCAAATGGG", Read(0, 10, "8M1S"))));
		g.addEvidence(SCE(FWD, withSequence("GCAAATGGG", Read(0, 10, "8M1S"))));
		g.addEvidence(SCE(FWD, withSequence("CCAAATGGG", Read(0, 10, "8M1S")))); // worst weight at the immediate next kmer
		g.addEvidence(SCE(FWD, withSequence("GGTACCCAAATGGG", Read(0, 10, "13M1S")))); // but best overall
		// should take the shorter assembly that has a longer soft clip
		PathGraphAssembler<DeBruijnSubgraphNode, DeBruijnPathNode<DeBruijnSubgraphNode>> pga = PGA(g, ap, "TGGG");
		List<List<DeBruijnSubgraphNode>> result = flatten(pga.assembleContigs());
		assertEquals(1, result.size());
		assertEquals("GGTACCCAAATGGG", S(g, result.get(0)));
	}
	@Test
	public void should_assemble_simple_sequence() {
		AssemblyParameters ap = new AssemblyParameters();
		ap.k = 4;
		ap.subgraphAssemblyTraversalMaximumBranchingFactor = 16;
		ap.maxContigsPerAssembly = 1000;
		ap.maxBaseMismatchForCollapse = 0;
		
		DeBruijnReadGraph g = RG(ap.k);
		g.addEvidence(SCE(FWD, withSequence("GGTACCCAAATGGG", Read(0, 10, "10M4S"))));
		PathGraphAssembler<DeBruijnSubgraphNode, DeBruijnPathNode<DeBruijnSubgraphNode>> pga = PGA(g, ap, "TGGG");
		List<List<DeBruijnSubgraphNode>> result = flatten(pga.assembleContigs());
		assertEquals(1, result.size());
		assertEquals("GGTACCCAAATGGG", S(g, result.get(0)));
	}
	/**
	 * Assembly back into reference is required for small-medium insertion detection
	 */
	@Test
	public void should_assemble_viable_breakpoint() {
		AssemblyParameters ap = new AssemblyParameters();
		ap.k = 3;
		ap.subgraphAssemblyTraversalMaximumBranchingFactor = 16;
		ap.maxContigsPerAssembly = 1000;
		ap.maxBaseMismatchForCollapse = 0;
		
		DeBruijnReadGraph g = RG(ap.k);
		g.addEvidence(SCE(FWD, withSequence("ACTGT", Read(0, 10, "3M2S"))));
		g.addEvidence(SCE(FWD, withSequence("GTCA", Read(0, 10, "3M1S"))));
		
		PathGraphAssembler<DeBruijnSubgraphNode, DeBruijnPathNode<DeBruijnSubgraphNode>> pga = PGA(g, ap, "ACT");
		List<List<DeBruijnSubgraphNode>> result = flatten(pga.assembleContigs());
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
		ap.subgraphAssemblyTraversalMaximumBranchingFactor = 1;
		ap.maxContigsPerAssembly = 4;
		ap.maxBaseMismatchForCollapse = 0;
		
		DeBruijnReadGraph g = RG(ap.k);
		//                                    RRRRRR
		g.addEvidence(SCE(FWD, withSequence( "TTTTCGAGAT", Read(0, 10, "8M1S"))));
		g.addEvidence(SCE(FWD, withSequence( "TTTTCGCCGGT", Read(0, 10, "8M3S"))));
		g.addEvidence(SCE(FWD, withSequence( "TTTTCGTTAACC", Read(0, 10, "8M6S"))));
		PathGraphAssembler<DeBruijnSubgraphNode, DeBruijnPathNode<DeBruijnSubgraphNode>> pga = PGA(g, ap,"TTTT");
		List<List<DeBruijnSubgraphNode>> result = flatten(pga.assembleContigs());
		assertEquals(3, result.size());
		assertEquals("TTTTCGTTAACC", S(g, result.get(0)));
		assertEquals("TTTTCGCCGGT", S(g, result.get(1)));
		assertEquals("TTTTCGAGAT", S(g, result.get(2)));
	}
	@Test
	public void should_fall_back_to_greedy_assemble() {
		AssemblyParameters ap = new AssemblyParameters();
		ap.k = 4;
		ap.maxContigsPerAssembly = 1000;
		ap.maxBaseMismatchForCollapse = 0;
		ap.maxPathTraversalNodes = 2;
		
		DeBruijnReadGraph g = RG(ap.k);
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
		PathGraphAssembler<DeBruijnSubgraphNode, DeBruijnPathNode<DeBruijnSubgraphNode>> pga = PGA(g, ap, "GGGG");
		//new StaticDeBruijnSubgraphPathGraphGexfExporter(ap.k)
		//	.snapshot(pga)
		//	.saveTo(new File("C:\\temp\\test.gexf"));
		List<List<DeBruijnSubgraphNode>> result = flatten(pga.assembleContigs());
		assertEquals(1, result.size());
		assertEquals("AAATGGGGA", S(g, result.get(0)));
	}
	@Test
	public void should_limit_anchor_assembly_length() {
		AssemblyParameters ap = new AssemblyParameters();
		ap.k = 5;
		ap.maxContigsPerAssembly = 1000;
		ap.maxBaseMismatchForCollapse = 0;
		ap.maxPathTraversalNodes = 2;
		ap.anchorAssemblyLength = 2;
		DeBruijnReadGraph g = RG(ap.k);
		g.addEvidence(SCE(FWD, withSequence( "AACCGGACAGCGCAGAGCATT", Read(0, 1, "20M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAGCGCAGAgGCATT", Read(0, 1, "20M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAGCGCAGgAGCATT", Read(0, 1, "20M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAGCGCAgGAGCATT", Read(0, 1, "20M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAGCGCgAGAGCATT", Read(0, 1, "20M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAGCGgCAGAGCATT", Read(0, 1, "20M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAGCgGCAGAGCATT", Read(0, 1, "20M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAGgCGCAGAGCATT", Read(0, 1, "20M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAgGCGCAGAGCATT", Read(0, 1, "20M1S"))));
		PathGraphAssembler<DeBruijnSubgraphNode, DeBruijnPathNode<DeBruijnSubgraphNode>> pga = PGA(g, ap, "GCATT");
		List<List<DeBruijnSubgraphNode>> result = flatten(pga.assembleContigs());
		assertEquals(1, result.size());
		assertEquals("GAGCATT", S(g, result.get(0))); 
		
		// above limit should fully assemble
		ap.anchorAssemblyLength = 1000;
		g = RG(ap.k);
		g.addEvidence(SCE(FWD, withSequence( "AACCGGACAGCGCAGAGCATT", Read(0, 1, "20M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAGCGCAGAgGCATT", Read(0, 1, "20M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAGCGCAGgAGCATT", Read(0, 1, "20M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAGCGCAgGAGCATT", Read(0, 1, "20M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAGCGCgAGAGCATT", Read(0, 1, "20M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAGCGgCAGAGCATT", Read(0, 1, "20M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAGCgGCAGAGCATT", Read(0, 1, "20M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAGgCGCAGAGCATT", Read(0, 1, "20M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAgGCGCAGAGCATT", Read(0, 1, "20M1S"))));
		pga = PGA(g, ap, "GCATT");
		result = flatten(pga.assembleContigs());
		assertEquals(1, result.size());
		assertEquals(17, result.get(0).size());
	}
}