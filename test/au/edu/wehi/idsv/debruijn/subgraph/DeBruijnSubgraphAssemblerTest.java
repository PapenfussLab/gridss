package au.edu.wehi.idsv.debruijn.subgraph;

import static org.junit.Assert.assertEquals;

import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.AssemblyParameters;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.VariantContextDirectedEvidence;

import com.google.common.collect.Lists;


public class DeBruijnSubgraphAssemblerTest extends TestHelper {
	private DeBruijnSubgraphAssembler DSA(int k) {
		ProcessingContext pc = getContext();
		AssemblyParameters p = pc.getAssemblyParameters();
		p.k = k;
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
}
