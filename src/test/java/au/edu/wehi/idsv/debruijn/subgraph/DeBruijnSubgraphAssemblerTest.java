package au.edu.wehi.idsv.debruijn.subgraph;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.util.List;

import org.junit.Ignore;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import au.edu.wehi.idsv.AssemblyParameters;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.VariantContextDirectedEvidence;

import com.google.common.collect.Lists;


public class DeBruijnSubgraphAssemblerTest extends TestHelper {
	@Rule
    public TemporaryFolder testFolder = new TemporaryFolder();
	private DeBruijnSubgraphAssembler DSA(int k) {
		ProcessingContext pc = getContext();
		AssemblyParameters p = pc.getAssemblyParameters();
		p.k = k;
		p.debruijnGraphVisualisationDirectory = new File(testFolder.getRoot(), "visualisation");
		p.visualiseAll = true;
		return new DeBruijnSubgraphAssembler(pc, AES());
	}
	@Test
	public void should_assemble_all_contigs() {
		DeBruijnSubgraphAssembler ass = DSA(3);
		List<VariantContextDirectedEvidence> results = Lists.newArrayList();
		results.addAll(Lists.newArrayList(ass.addEvidence(SCE(FWD, withSequence("TAAAGTC", Read(0, 1, "4M3S"))))));
		results.addAll(Lists.newArrayList(ass.addEvidence(SCE(FWD, withSequence("AAAGTCT", Read(0, 2, "3M4S"))))));
		results.addAll(Lists.newArrayList(ass.addEvidence(SCE(FWD, withSequence("TAAAGTC", Read(1, 1, "4M3S"))))));
		results.addAll(Lists.newArrayList(ass.addEvidence(SCE(FWD, withSequence("AAAGTCT", Read(1, 2, "3M4S"))))));
		results.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(2, results.size());
	}
	@Test
	public void should_export_debruijn_graph() {
		DeBruijnSubgraphAssembler ass = DSA(3);
		List<VariantContextDirectedEvidence> results = Lists.newArrayList();
		results.addAll(Lists.newArrayList(ass.addEvidence(SCE(FWD, withSequence("TAAAGTC", Read(0, 1, "4M3S"))))));
		results.addAll(Lists.newArrayList(ass.addEvidence(SCE(FWD, withSequence("AAAGTCT", Read(0, 2, "3M4S"))))));
		results.addAll(Lists.newArrayList(ass.addEvidence(SCE(FWD, withSequence("TAAAGTC", Read(1, 1, "4M3S"))))));
		results.addAll(Lists.newArrayList(ass.addEvidence(SCE(FWD, withSequence("AAAGTCT", Read(1, 2, "3M4S"))))));
		results.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertTrue(new File(new File(testFolder.getRoot(), "visualisation"), "debruijn.kmers.forward.polyA.gexf").exists());
		assertTrue(new File(new File(testFolder.getRoot(), "visualisation"), "debruijn.kmers.backward.polyA.gexf").exists());
		assertTrue(new File(new File(testFolder.getRoot(), "visualisation"), "debruijn.kmers.forward.polyACGT.gexf").exists());
		assertTrue(new File(new File(testFolder.getRoot(), "visualisation"), "debruijn.kmers.backward.polyACGT.gexf").exists());
	}
	@Test
	public void should_assemble_both_directions() {
		DeBruijnSubgraphAssembler ass = DSA(3);
		List<VariantContextDirectedEvidence> results = Lists.newArrayList();
		results.addAll(Lists.newArrayList(ass.addEvidence(SCE(FWD, withSequence("TAAAGTC", Read(0, 1, "4M3S"))))));
		results.addAll(Lists.newArrayList(ass.addEvidence(SCE(FWD, withSequence("AAAGTCT", Read(0, 2, "3M4S"))))));
		results.addAll(Lists.newArrayList(ass.addEvidence(SCE(BWD, withSequence("TATG", Read(0, 10, "1S3M"))))));
		results.addAll(Lists.newArrayList(ass.addEvidence(SCE(BWD, withSequence("TTATG", Read(0, 10, "2S3M"))))));
		results.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(2, results.size());
	}
	@Test
	public void should_anchor_at_reference_kmer() {
		DeBruijnSubgraphAssembler ass = DSA(3);
		List<VariantContextDirectedEvidence> results = Lists.newArrayList();
		results.addAll(Lists.newArrayList(ass.addEvidence(SCE(FWD, withSequence("TAAAGTC", Read(0, 1, "4M3S"))))));
		results.addAll(Lists.newArrayList(ass.addEvidence(SCE(FWD, withSequence("AAAGTCT", Read(0, 2, "3M4S"))))));
		results.addAll(Lists.newArrayList(ass.addEvidence(SCE(BWD, withSequence("CTGAAAT", Read(0, 10, "3S4M"))))));
		results.addAll(Lists.newArrayList(ass.addEvidence(SCE(BWD, withSequence("TCTGAAA", Read(0, 10, "4S3M"))))));
		results.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(2, results.size());
		assertEquals(4, results.get(0).getBreakendSequence().length);
	}
	@Test
	public void should_anchor_at_reference_kmer_large_kmer() {
		DeBruijnSubgraphAssembler ass = DSA(32);
		List<VariantContextDirectedEvidence> results = Lists.newArrayList();
		results.addAll(Lists.newArrayList(ass.addEvidence(SCE(FWD, withSequence(S(RANDOM).substring(0, 200), Read(0, 1, "100M100S"))))));
		results.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(100, results.get(0).getBreakendSequence().length);
	}
	@Test
	public void soft_clip_assembly_should_anchor_at_reference_kmer() {
		DeBruijnSubgraphAssembler ass = DSA(4);
		List<VariantContextDirectedEvidence> results = Lists.newArrayList();
		results.addAll(Lists.newArrayList(ass.addEvidence(SCE(BWD, withSequence("TTGCTCAAAA", Read(0, 1, "6S4M"))))));
		results.addAll(Lists.newArrayList(ass.addEvidence(NRRP(withSequence("TGCTG", OEA(0, 4, "5M", false))))));
		results.addAll(Lists.newArrayList(ass.addEvidence(NRRP(withSequence("TGCTG", OEA(0, 4, "5M", false))))));
		results.addAll(Lists.newArrayList(ass.addEvidence(NRRP(withSequence("TGCTG", OEA(0, 4, "5M", false))))));
		results.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(1, results.size());
		assertEquals("TTGCTCAAAA", results.get(0).getAssemblyConsensus());
		assertEquals(1, results.get(0).getBreakendSummary().start);
		assertEquals(4, results.get(0).getAssemblyConsensus().length() - results.get(0).getBreakendSequence().length);
	}
	@Test
	@Ignore("TODO: NYI: Not Yet Implemented")
	public void should_assemble_anchor_shorter_than_kmer() {
		DeBruijnSubgraphAssembler ass = DSA(5);
		List<VariantContextDirectedEvidence> results = Lists.newArrayList();
		results.addAll(Lists.newArrayList(ass.addEvidence(SCE(FWD, withSequence("ATTAGA", Read(0, 1, "1M5S"))))));
		results.addAll(Lists.newArrayList(ass.addEvidence(SCE(FWD, withSequence("ATTAGA", Read(0, 1, "1M5S"))))));
		results.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(1, results.size());
	}
	@Test
	@Ignore("TODO: NYI: Not Yet Implemented")
	public void should_assemble_anchor_shorter_than_kmer_with_indel_rp_support() {
		DeBruijnSubgraphAssembler ass = DSA(5);
		List<VariantContextDirectedEvidence> results = Lists.newArrayList();
		results.addAll(Lists.newArrayList(ass.addEvidence(NRRP(withSequence("GTCTTA", DP(0, 1, "8M", true, 0, 500, "8M", false))))));
		results.addAll(Lists.newArrayList(ass.addEvidence(SCE(FWD, withSequence("CTTAGA", Read(0, 100, "1M5S"))))));
		results.addAll(Lists.newArrayList(ass.addEvidence(SCE(FWD, withSequence("CTTAGA", Read(0, 100, "1M5S"))))));
		results.addAll(Lists.newArrayList(ass.endOfEvidence()));
		assertEquals(1, results.size());
		assertEquals(3, results.get(0).getBreakendSummary().start);
	}
}
