package au.edu.wehi.idsv.debruijn.subgraph;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.AssemblyParameters;
import au.edu.wehi.idsv.AssemblyParameters.ContigAssemblyOrder;
import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.VariantContextDirectedEvidence;
import au.edu.wehi.idsv.vcf.VcfAttributes;
import au.edu.wehi.idsv.vcf.VcfFilter;

import com.google.common.collect.Lists;


public class DeBruijnReadGraphTest extends TestHelper {
	public DeBruijnReadGraph G(int referenceIndex, int k, BreakendDirection direction) {
		AssemblyParameters p = new AssemblyParameters();
		p.maxBaseMismatchForCollapse = 0;
		p.assemblyOrder = ContigAssemblyOrder.GreedyMaxKmer;
		p.k = k;
		return new DeBruijnReadGraph(getContext(), AES(), referenceIndex, direction, p);
	}
	@Test
	public void should_assemble_sc() {
		DeBruijnReadGraph g = G(0, 3, FWD);
		g.addEvidence(SCE(FWD, withSequence("TAAAGTC", Read(0, 10, "4M3S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAGTCT", Read(0, 11, "3M4S"))));
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		assertEquals(1, result.size());
		VariantContextDirectedEvidence bp = result.get(0);
		assertEquals("TAAAGTCT", bp.getAssemblyConsensus());
	}
	@Test
	public void should_assemble_sc_b() {
		DeBruijnReadGraph g = G(0, 3, BWD);
		g.addEvidence(SCE(BWD, withSequence("TATG", Read(0, 10, "1S3M"))));
		g.addEvidence(SCE(BWD, withSequence("TTATG", Read(0, 10, "2S3M"))));
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		assertEquals(1, result.size());
		VariantContextDirectedEvidence bp = result.get(0);
		assertEquals("TTATG", bp.getAssemblyConsensus());
	}
	@Test
	public void should_set_assembly_common_attributes() {
		DeBruijnReadGraph g = G(0, 3, FWD);
		g.addEvidence(SCE(FWD, withQual(new byte[] { 4,4,4,4,4,4,4}, withSequence("AAAGTCT", Read(0, 10, "3M4S")))));
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		VariantContextDirectedEvidence bp = result.get(0);
		assertEquals("AAAGTCT", bp.getAssemblyConsensus());
		assertEquals(4, bp.getAssemblyMaximumSoftClipLength());
		assertEquals("debruijn-s", bp.getAssemblerProgram());
		assertEquals(4d, bp.getAssemblyQuality(), 0);
	}
	@Test
	public void should_set_assembly_sc_attributes() {
		DeBruijnReadGraph g = G(0, 3, FWD);
		g.addEvidence(SCE(FWD, withQual(new byte[] { 4,4,4,4,4,4,4}, withSequence("TAAAGTC", Read(0, 10, "4M3S")))));
		g.addEvidence(SCE(FWD, withQual(new byte[] { 4,4,4,4,4,4,4}, withSequence("AAAGTCT", Read(0, 10, "3M4S")))));
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		VariantContextDirectedEvidence bp = result.get(0);
		assertEquals(4, bp.getAssemblyMaximumSoftClipLength());
	}
	public void should_set_assembly_mate_attributes() {
		DeBruijnReadGraph g = G(0, 3, FWD);
		g.addEvidence(NRRP(OEA(0, 1, "5M", true)));
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		VariantContextDirectedEvidence bp = result.get(0);
		assertEquals(5, bp.getAssemblyLongestSupportingRead());
	}
	@Test
	public void should_set_evidence_attributes() {
		DeBruijnReadGraph g = G(0, 3, FWD);
		g.addEvidence(SCE(FWD, withSequence("TAAAGTC", Read(0, 10, "4M3S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAGTCT", Read(0, 10, "3M4S"))));
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		VariantContextDirectedEvidence bp = result.get(0);
		assertEquals(bp.getAssemblyConsensus().length(), bp.getBreakendSummary().evidence.get(VcfAttributes.ASSEMBLY_LENGTH));
		assertEquals(2, bp.getBreakendSummary().evidence.get(VcfAttributes.ASSEMBLY_READS));
		assertEquals(14, bp.getBreakendSummary().evidence.get(VcfAttributes.ASSEMBLY_BASES));
	}
	@Test
	public void should_assemble_adjacent_scs() {
		DeBruijnReadGraph g = G(0, 3, FWD);
		g.addEvidence(SCE(FWD, withSequence("AAAGTC", Read(0, 10, "3M3S")))); // 12
		g.addEvidence(SCE(FWD, withSequence("AAAAGTCTT", Read(0, 10, "4M4S")))); // 13
		g.addEvidence(SCE(FWD, withSequence("AAAAGTCTT", Read(0, 10, "4M4S")))); // 13
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		assertEquals(1, result.size());
		VariantContextDirectedEvidence bp = result.get(0);
		assertEquals("AAAGTCTT", bp.getAssemblyConsensus());
		assertEquals(0, bp.getBreakendSummary().referenceIndex);
		assertEquals(13, bp.getBreakendSummary().start);
		assertEquals(13, bp.getBreakendSummary().end);
	}
	@Test
	public void should_anchor_to_best_reference_position() {
		DeBruijnReadGraph g = G(0, 3, FWD);
		g.addEvidence(SCE(FWD, withSequence("AAAGTC", Read(0, 9, "3M3S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAGTC", Read(0, 10, "3M3S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAGTC", Read(0, 10, "3M3S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAGTC", Read(0, 10, "3M3S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAGTC", Read(0, 11, "3M3S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAAGTCT", Read(0, 10, "4M4S"))));
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		// we should call to 12 since we have many reads supporting that assembly
		assertEquals(12, result.get(0).getBreakendSummary().start);
	}
	@Test
	public void should_assemble_sc_and_nrrp_together() {
		DeBruijnReadGraph g = G(0, 4, FWD);
		g.addEvidence(SCE(FWD, withSequence("TAAAGTCC", Read(0, 10, "4M4S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAAGTCCT", Read(0, 10, "4M5S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAAGTCCT", Read(0, 10, "4M5S"))));
		g.addEvidence(NRRP(withSequence("GTCCTAGAC", DP(0, 1, "8M", true, 1, 10, "8M", false))));
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		assertEquals(1, result.size());
		VariantContextDirectedEvidence bp = result.get(0);
		//       GTCC <-seed
		//        TCCT
		//         CCTA
		//          CTAG
		//           TAGA
		//            AGAC
		//      AGTC
		//     AAGT 
		//    AAAG 
		//   AAAA
		//  TAAA  
		assertEquals("TAAAAGTCCTAGAC", bp.getAssemblyConsensus());
	}
	@Test
	public void should_assemble_nrrp() {
		DeBruijnReadGraph g = G(0, 3, FWD);
		// expect FR orientation so we need to reverse comp the unmapped seq
		g.addEvidence(NRRP(withSequence("AGAC", OEA(0, 10, "4M", true)))); // revcomp=GTCT
		g.addEvidence(NRRP(withSequence("CTAG", DP(0, 9, "4M", true, 1, 10, "4M", false))));
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		assertEquals(1, result.size());
		VariantContextDirectedEvidence bp = result.get(0);
		assertEquals("GTCTAG", bp.getAssemblyConsensus());
		assertEquals(FWD, bp.getBreakendSummary().direction);
		assertEquals(0, bp.getBreakendSummary().referenceIndex);
		assertEquals(13, bp.getBreakendSummary().start);
		assertEquals(313, bp.getBreakendSummary().end);
	}
	@Test
	public void should_assemble_best_contig() {
		DeBruijnReadGraph g = G(0, 3, FWD);
		g.addEvidence(SCE(FWD, withSequence("AAAGTCTA", Read(0, 10, "3M5S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAGTCTA", Read(0, 10, "3M5S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAGTCTG", Read(0, 10, "3M5S"))));
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		assertEquals(1, result.size());
		VariantContextDirectedEvidence bp = result.get(0);
		// CTG becomes unanchored & we can't make a contig out of it
		assertEquals("AAAGTCTA", bp.getAssemblyConsensus());
	}
	@Test
	public void should_assemble_when_out_of_scope() {
		DeBruijnReadGraph g = G(0, 3, FWD);
		g.addEvidence(SCE(FWD, withSequence("TAAAGTC", Read(0, 1, "4M3S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAGTCT", Read(0, 2, "3M4S"))));
		assertEquals(0, Lists.newArrayList(g.assembleContigsBefore(3)).size());
		assertEquals(0, Lists.newArrayList(g.assembleContigsBefore(4)).size()); // anchored at position 4
		assertEquals(1, Lists.newArrayList(g.assembleContigsBefore(5)).size()); // so should subgraph after this position
	}
	@Test
	public void should_filter_short_contigs_that_do_not_include_reference_bases() {
		DeBruijnReadGraph g = G(0, 4, FWD);
		g.addEvidence(SCE(FWD, withSequence("TAAAGTCC", Read(0, 10, "4M4S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAAGTCCT", Read(0, 10, "4M5S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAAGTCCT", Read(0, 10, "4M5S"))));
		g.addEvidence(NRRP(withSequence("GTCCTAGAC", DP(0, 1, "8M", true, 1, 10, "8M", false))));
		g.addEvidence(NRRP(withSequence("GTCCTAGAC", DP(0, 1, "8M", true, 1, 10, "8M", false))));
		g.addEvidence(NRRP(withSequence("GTCCTAGAC", DP(0, 1, "8M", true, 1, 10, "8M", false))));
		g.addEvidence(NRRP(withSequence("GTCCTAGAT", DP(0, 1, "8M", true, 1, 10, "8M", false))));
		g.addEvidence(NRRP(withSequence("GTCCTAGAT", DP(0, 1, "8M", true, 1, 10, "8M", false))));
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		assertEquals(2, result.size());
		assertEquals("TAAAAGTCCTAGAC", result.get(1).getAssemblyConsensus());
		// AGAT is left over as a seed but since it does not anchor to the reference
		// As it's shorter than a read assembly, we presume that it is
		// sequencing noise and abandon it
		assertEquals("AGAT", result.get(0).getAssemblyConsensus());
		assertTrue(result.get(0).getFilters().contains(VcfFilter.ASSEMBLY_TOO_SHORT.name()));
	}
	@Test
	public void debugPrintPaths_should_work() {
		DeBruijnReadGraph g = G(0, 4, FWD);
		g.addEvidence(SCE(FWD, withSequence("TAAAGTCC", Read(0, 10, "4M4S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAAGTCCT", Read(0, 10, "4M5S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAAGTCCT", Read(0, 10, "4M5S"))));
		g.addEvidence(NRRP(withSequence("GTCCTAGAC", DP(0, 1, "8M", true, 1, 10, "8M", false))));
		g.addEvidence(NRRP(withSequence("GTCCTAGAC", DP(0, 1, "8M", true, 1, 10, "8M", false))));
		g.addEvidence(NRRP(withSequence("GTCCTAGAC", DP(0, 1, "8M", true, 1, 10, "8M", false))));
		g.addEvidence(NRRP(withSequence("GTCCTAGAT", DP(0, 1, "8M", true, 1, 10, "8M", false))));
		g.addEvidence(NRRP(withSequence("GTCCTAGAT", DP(0, 1, "8M", true, 1, 10, "8M", false))));
		String s = g.debugPrintPaths();
		assertNotNull(s);
	}
	@Test
	public void reference_kmer_support_should_not_be_considered_read_support() {
		DeBruijnReadGraph g = G(0, 3, FWD);
		g.addEvidence(SCE(FWD, withSequence("AAAAAAAA", Read(0, 10, "3M5S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAGTCTA", Read(0, 10, "3M5S"))));
		List<VariantContextDirectedEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		assertEquals(1, result.size());
		VariantContextDirectedEvidence bp = result.get(0);
		// CTG becomes unanchored & we can't make a contig out of it
		assertEquals(1, bp.getBreakendSummary().evidence.get(VcfAttributes.ASSEMBLY_READS));
	}
}
