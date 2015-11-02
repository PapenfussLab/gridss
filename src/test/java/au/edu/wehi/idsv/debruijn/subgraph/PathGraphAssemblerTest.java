package au.edu.wehi.idsv.debruijn.subgraph;

import static org.junit.Assert.assertEquals;

import java.util.List;

import org.junit.Ignore;
import org.junit.Test;

import com.google.common.base.Function;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;

import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.configuration.AssemblyConfiguration;
import au.edu.wehi.idsv.configuration.GridssConfiguration;
import au.edu.wehi.idsv.debruijn.DeBruijnPathNode;
import au.edu.wehi.idsv.debruijn.DeBruijnPathNodeFactory;
import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.graph.PathNode;
import au.edu.wehi.idsv.visualisation.NontrackingSubgraphTracker;


public class PathGraphAssemblerTest extends TestHelper {
	PathGraphAssembler<DeBruijnSubgraphNode, DeBruijnPathNode<DeBruijnSubgraphNode>> PGA(DeBruijnReadGraph g, GridssConfiguration config, String seed) {
		DeBruijnSubgraphNode seedNode = g.getKmer(KmerEncodingHelper.picardBaseToEncoded(config.getAssembly().k, B(seed)));
		PathGraphAssembler<DeBruijnSubgraphNode, DeBruijnPathNode<DeBruijnSubgraphNode>> pga = new PathGraphAssembler<DeBruijnSubgraphNode, DeBruijnPathNode<DeBruijnSubgraphNode>>(
			g,
			config,
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
		GridssConfiguration config = getConfig();
		AssemblyConfiguration ap = config.getAssembly();
		ap.k = 4;
		ap.subgraph.traversalMaximumBranchingFactor = 1;
		ap.subgraph.maxContigsPerAssembly = 1;
		
		DeBruijnReadGraph g = RG(ap.k);
		g.addEvidence(SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))));
		g.addEvidence(SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))));
		g.addEvidence(NRRP(withSequence("GGTTAACC", DP(0, 1, "8M", true, 1, 10, "8M", false))));
		g.addEvidence(SCE(FWD, withSequence("TTGGT", Read(0, 10, "4M1S"))));
		//     TTAACCGGCCAATT
		//   GGTTAACC
		// TTGGT
		// ^   ^^^^ <- starts of reference kmers
		PathGraphAssembler<DeBruijnSubgraphNode, DeBruijnPathNode<DeBruijnSubgraphNode>> pga = PGA(g, config, "AATT");
		List<List<DeBruijnSubgraphNode>> result = flatten(pga.assembleContigs());
		assertEquals(2, result.size());
		assertEquals("TTAACCGGCCAATT", S(g, result.get(0)));
	}
	@Ignore("now removing entire non-reference subgraph after each assembly")
	//@Test
	public void should_assemble_excluding_used_non_reference_kmers() {
		GridssConfiguration config = getConfig();
		AssemblyConfiguration ap = config.getAssembly();
		ap.k = 4;
		ap.subgraph.traversalMaximumBranchingFactor = 1;
		ap.subgraph.maxContigsPerAssembly = 1000;
		ap.errorCorrection.maxBaseMismatchForCollapse = 0;
		
		DeBruijnReadGraph g = RG(ap.k);
		g.addEvidence(SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))));
		g.addEvidence(SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))));
		g.addEvidence(NRRP(withSequence(           "GCCAATC", OEA(0, 10, "7M", true)))); // not on main patch, but reachable
		g.addEvidence(SCE(FWD, withSequence("TTAACCGAGTCCTG", Read(0, 10, "7M7S"))));
		
		PathGraphAssembler<DeBruijnSubgraphNode, DeBruijnPathNode<DeBruijnSubgraphNode>> pga = PGA(g, config, "AATT");
		List<List<DeBruijnSubgraphNode>> result = flatten(pga.assembleContigs());
		assertEquals(3, result.size());
		assertEquals("TTAACCGGCCAATT", S(g, result.get(0)));
		assertEquals("TTAACCGAGTCCTG", S(g, result.get(1)));
	}
	@Test
	public void should_assemble_excluding_all_reachable_non_reference_kmers() {
		GridssConfiguration config = getConfig();
		AssemblyConfiguration ap = config.getAssembly();
		ap.k = 4;
		ap.subgraph.traversalMaximumBranchingFactor = 1;
		ap.subgraph.maxContigsPerAssembly = 1000;
		ap.errorCorrection.maxBaseMismatchForCollapse = 0;
		
		
		DeBruijnReadGraph g = RG(ap.k);
		g.addEvidence(SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))));
		g.addEvidence(SCE(FWD, withSequence("TTAACCGGCCAATT", Read(0, 10, "7M7S"))));
		g.addEvidence(NRRP(withSequence(              "AATC", OEA(0, 10, "4M", true)))); // not on main path, but reachable
		g.addEvidence(SCE(FWD, withSequence("TTAACCGAGTCCTG", Read(0, 10, "7M7S"))));
		
		PathGraphAssembler<DeBruijnSubgraphNode, DeBruijnPathNode<DeBruijnSubgraphNode>> pga = PGA(g, config, "AATT");
		List<List<DeBruijnSubgraphNode>> result = flatten(pga.assembleContigs());
		assertEquals(2, result.size());
		assertEquals("TTAACCGGCCAATT", S(g, result.get(0)));
		assertEquals("TTAACCGAGTCCTG", S(g, result.get(1)));
	}
	@Test
	public void should_maximise_non_reference_weight() {
		GridssConfiguration config = getConfig();
		AssemblyConfiguration ap = config.getAssembly();
		ap.k = 4;
		ap.subgraph.maxContigsPerAssembly = 1000;
		ap.errorCorrection.maxBaseMismatchForCollapse = 0;
		
		DeBruijnReadGraph g = RG(ap.k);
		g.addEvidence(SCE(FWD, withSequence("AAATGGGG", Read(0, 10, "6M2S"))));
		g.addEvidence(SCE(FWD, withSequence("GTACCCGGGG", Read(0, 10, "9M1S"))));
		// should take the shorter assembly that has a longer soft clip
		PathGraphAssembler<DeBruijnSubgraphNode, DeBruijnPathNode<DeBruijnSubgraphNode>> pga = PGA(g, config, "GGGG");
		List<List<DeBruijnSubgraphNode>> result = flatten(pga.assembleContigs());
		assertEquals(1, result.size());
		assertEquals("AAATGGGG", S(g, result.get(0)));
	}
	@Test
	public void greedy_should_maximise_reference_weight() {
		GridssConfiguration config = getConfig();
		AssemblyConfiguration ap = config.getAssembly();
		ap.k = 4;
		ap.subgraph.traversalMaximumBranchingFactor = 1;
		ap.subgraph.maxContigsPerAssembly = 1000;
		ap.errorCorrection.maxBaseMismatchForCollapse = 0;
		
		DeBruijnReadGraph g = RG(ap.k);
		g.addEvidence(SCE(FWD, withSequence( "TAAATGGG", Read(0, 10, "7M1S"))));
		g.addEvidence(SCE(FWD, withSequence( "TAAATGGG", Read(0, 10, "7M1S"))));
		g.addEvidence(SCE(FWD, withSequence( "TAAATGGG", Read(0, 10, "7M1S"))));
		g.addEvidence(SCE(FWD, withSequence("GCAAATGGG", Read(0, 10, "8M1S"))));
		g.addEvidence(SCE(FWD, withSequence("GCAAATGGG", Read(0, 10, "8M1S"))));
		g.addEvidence(SCE(FWD, withSequence("GCAAATGGG", Read(0, 10, "8M1S"))));
		g.addEvidence(SCE(FWD, withSequence("CCAAATGGG", Read(0, 10, "8M1S")))); // worst weight at the immediate next kmer
		g.addEvidence(SCE(FWD, withSequence("GGTACCCAAATGGG", Read(0, 10, "13M1S")))); // but best overall path node
		// should take the shorter assembly that has a longer soft clip
		PathGraphAssembler<DeBruijnSubgraphNode, DeBruijnPathNode<DeBruijnSubgraphNode>> pga = PGA(g, config, "TGGG");
		List<List<DeBruijnSubgraphNode>> result = flatten(pga.assembleContigs());
		assertEquals(1, result.size());
		assertEquals("GGTACCCAAATGGG", S(g, result.get(0)));
	}
	@Test
	public void should_assemble_simple_sequence() {
		GridssConfiguration config = getConfig();
		AssemblyConfiguration ap = config.getAssembly();
		ap.k = 4;
		ap.subgraph.traversalMaximumBranchingFactor = 16;
		ap.subgraph.maxContigsPerAssembly = 1000;
		ap.errorCorrection.maxBaseMismatchForCollapse = 0;
		
		DeBruijnReadGraph g = RG(ap.k);
		g.addEvidence(SCE(FWD, withSequence("GGTACCCAAATGGG", Read(0, 10, "10M4S"))));
		PathGraphAssembler<DeBruijnSubgraphNode, DeBruijnPathNode<DeBruijnSubgraphNode>> pga = PGA(g, config, "TGGG");
		List<List<DeBruijnSubgraphNode>> result = flatten(pga.assembleContigs());
		assertEquals(1, result.size());
		assertEquals("GGTACCCAAATGGG", S(g, result.get(0)));
	}
	/**
	 * Assembly back into reference is required for small-medium insertion detection
	 */
	@Test
	public void should_assemble_viable_breakpoint() {
		GridssConfiguration config = getConfig();
		AssemblyConfiguration ap = config.getAssembly();
		ap.k = 3;
		ap.subgraph.traversalMaximumBranchingFactor = 16;
		ap.subgraph.maxContigsPerAssembly = 1000;
		ap.errorCorrection.maxBaseMismatchForCollapse = 0;
		
		DeBruijnReadGraph g = RG(ap.k);
		g.addEvidence(SCE(FWD, withSequence("ACTGT", Read(0, 10, "3M2S"))));
		g.addEvidence(SCE(FWD, withSequence("GTCA", Read(0, 10, "3M1S"))));
		
		PathGraphAssembler<DeBruijnSubgraphNode, DeBruijnPathNode<DeBruijnSubgraphNode>> pga = PGA(g, config, "ACT");
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
		GridssConfiguration config = getConfig();
		AssemblyConfiguration ap = config.getAssembly();
		ap.k = 4;
		ap.subgraph.traversalMaximumBranchingFactor = 1;
		ap.subgraph.maxContigsPerAssembly = 4;
		ap.errorCorrection.maxBaseMismatchForCollapse = 0;
		
		DeBruijnReadGraph g = RG(ap.k);
		//                                    RRRRRR
		g.addEvidence(SCE(FWD, withSequence( "TTTTCGAGAT", Read(0, 10, "8M2S"))));
		g.addEvidence(SCE(FWD, withSequence( "TTTTCGCCGGT", Read(0, 10, "8M3S"))));
		g.addEvidence(SCE(FWD, withSequence( "TTTTCGTTAACC", Read(0, 10, "8M4S"))));
		PathGraphAssembler<DeBruijnSubgraphNode, DeBruijnPathNode<DeBruijnSubgraphNode>> pga = PGA(g, config,"TTTT");
		List<List<DeBruijnSubgraphNode>> result = flatten(pga.assembleContigs());
		assertEquals(3, result.size());
		assertEquals("TTTTCGTTAACC", S(g, result.get(0)));
		assertEquals("TTTTCGCCGGT", S(g, result.get(1)));
		assertEquals("TTTTCGAGAT", S(g, result.get(2)));
	}
	@Test
	public void should_fall_back_to_greedy_assemble() {
		GridssConfiguration config = getConfig();
		AssemblyConfiguration ap = config.getAssembly();
		ap.k = 4;
		ap.subgraph.maxContigsPerAssembly = 1000;
		ap.errorCorrection.maxBaseMismatchForCollapse = 0;
		ap.subgraph.traveralMaximumPathNodes = 2;
		
		DeBruijnReadGraph g = RG(ap.k);
		g.addEvidence(SCE(FWD, withSequence("AAATGGGGGGGGGGG", Read(0, 10, "6M9S"))));
		g.addEvidence(SCE(FWD, withSequence("AAATGGGC", Read(0, 10, "6M2S"))));
		g.addEvidence(SCE(FWD, withSequence("AAATGGGTAA", Read(0, 10, "6M4S"))));
		g.addEvidence(SCE(FWD, withSequence("AAATGGGTCC", Read(0, 10, "6M4S"))));
		g.addEvidence(NRRP(withSequence("GGAC", OEA(0, 10, "4M", true)))); // GTCC reverse comp GGAC
		g.addEvidence(NRRP(withSequence("GGAC", OEA(0, 10, "4M", true))));
		g.addEvidence(SCE(FWD, withSequence("AAATGGGA", Read(0, 10, "6M2S"))));
		g.addEvidence(SCE(FWD, withSequence("AAATGGGA", Read(0, 10, "6M2S"))));
		g.addEvidence(SCE(FWD, withSequence("AAATGGGA", Read(0, 10, "6M2S"))));
		g.addEvidence(SCE(FWD, withSequence("AAATGGGA", Read(0, 10, "6M2S"))));
		// greedy traversal
		PathGraphAssembler<DeBruijnSubgraphNode, DeBruijnPathNode<DeBruijnSubgraphNode>> pga = PGA(g, config, "GGGG");
		//new StaticDeBruijnSubgraphPathGraphGexfExporter(ap.k)
		//	.snapshot(pga)
		//	.saveTo(new File("C:\\temp\\test.gexf"));
		List<List<DeBruijnSubgraphNode>> result = flatten(pga.assembleContigs());
		assertEquals(1, result.size());
		assertEquals("AAATGGGGA", S(g, result.get(0)));
	}
	@Test
	public void should_limit_anchor_assembly_length() {
		GridssConfiguration config = getConfig();
		AssemblyConfiguration ap = config.getAssembly();
		ap.k = 5;
		ap.subgraph.maxContigsPerAssembly = 1000;
		ap.errorCorrection.maxBaseMismatchForCollapse = 0;
		ap.subgraph.traveralMaximumPathNodes = 2;
		ap.anchorLength = 2;
		DeBruijnReadGraph g = RG(ap.k);
		g.addEvidence(SCE(FWD, withSequence( "AACCGGACAGCGCAGAGCATT", Read(0, 1, "20M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAGCGCAGAgGCATT", Read(0, 1, "21M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAGCGCAGgAGCATT", Read(0, 1, "21M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAGCGCAgGAGCATT", Read(0, 1, "21M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAGCGCgAGAGCATT", Read(0, 1, "21M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAGCGgCAGAGCATT", Read(0, 1, "21M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAGCgGCAGAGCATT", Read(0, 1, "21M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAGgCGCAGAGCATT", Read(0, 1, "21M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAgGCGCAGAGCATT", Read(0, 1, "21M1S"))));
		PathGraphAssembler<DeBruijnSubgraphNode, DeBruijnPathNode<DeBruijnSubgraphNode>> pga = PGA(g, config, "GCATT");
		List<List<DeBruijnSubgraphNode>> result = flatten(pga.assembleContigs());
		assertEquals(1, result.size());
		assertEquals("GAGCATT", S(g, result.get(0))); 
		
		// above limit should fully assemble
		ap.anchorLength = 1000;
		g = RG(ap.k);
		g.addEvidence(SCE(FWD, withSequence( "AACCGGACAGCGCAGAGCATT", Read(0, 1, "20M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAGCGCAGAgGCATT", Read(0, 1, "21M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAGCGCAGgAGCATT", Read(0, 1, "21M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAGCGCAgGAGCATT", Read(0, 1, "21M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAGCGCgAGAGCATT", Read(0, 1, "21M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAGCGgCAGAGCATT", Read(0, 1, "21M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAGCgGCAGAGCATT", Read(0, 1, "21M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAGgCGCAGAGCATT", Read(0, 1, "21M1S"))));
		g.addEvidence(SCE(FWD, withSequence("AACCGGACAgGCGCAGAGCATT", Read(0, 1, "21M1S"))));
		pga = PGA(g, config, "GCATT");
		result = flatten(pga.assembleContigs());
		assertEquals(1, result.size());
		assertEquals(17, result.get(0).size());
	}
}