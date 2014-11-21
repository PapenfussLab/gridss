package au.edu.wehi.idsv.debruijn.subgraph;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.AssemblyEvidence;
import au.edu.wehi.idsv.AssemblyParameters;
import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.EvidenceSubset;
import au.edu.wehi.idsv.TestHelper;

import com.google.common.collect.Lists;


public class DeBruijnReadGraphTest extends TestHelper {
	public DeBruijnReadGraph G(int referenceIndex, int k, BreakendDirection direction) {
		AssemblyParameters p = new AssemblyParameters();
		p.maxBaseMismatchForCollapse = 0;
		p.subgraphAssemblyTraversalMaximumBranchingFactor = 1;
		p.k = k;
		return new DeBruijnReadGraph(getContext(), AES(), referenceIndex, direction, p, null);
	}
	@Test
	public void should_excluse_kmers_containing_ambiguous_base() {
		DeBruijnReadGraph g = G(0, 3, FWD);
		g.addEvidence(SCE(FWD, withSequence("TANAGTN", Read(0, 10, "4M3S"))));
		g.addEvidence(SCE(FWD, withSequence( "AANGTCT", Read(0, 11, "3M4S"))));
		List<AssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		assertEquals(0, result.size());
	}
	@Test
	public void should_assemble_sc() {
		DeBruijnReadGraph g = G(0, 3, FWD);
		g.addEvidence(SCE(FWD, withSequence("TAAAGTC", Read(0, 10, "4M3S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAGTCT", Read(0, 11, "3M4S"))));
		List<AssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		assertEquals(1, result.size());
		AssemblyEvidence bp = result.get(0);
		assertEquals("TAAAGTCT", S(bp.getAssemblySequence()));
	}
	@Test
	public void should_assemble_sc_b() {
		DeBruijnReadGraph g = G(0, 3, BWD);
		g.addEvidence(SCE(BWD, withSequence("TATG", Read(0, 10, "1S3M"))));
		g.addEvidence(SCE(BWD, withSequence("TTATG", Read(0, 10, "2S3M"))));
		List<AssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		assertEquals(1, result.size());
		AssemblyEvidence bp = result.get(0);
		assertEquals("TTATG", S(bp.getAssemblySequence()));
	}
	@Test
	public void should_set_assembly_common_attributes() {
		DeBruijnReadGraph g = G(0, 3, FWD);
		g.addEvidence(SCE(FWD, withQual(new byte[] { 4,4,4,4,4,4,4}, withSequence("AAAGTCT", Read(0, 10, "3M4S")))));
		List<AssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		AssemblyEvidence bp = result.get(0);
		assertEquals("AAAGTCT", S(bp.getAssemblySequence()));
		assertEquals(4, bp.getAssemblySoftClipLengthMax(null));
	}
	@Test
	public void should_set_assembly_sc_attributes() {
		DeBruijnReadGraph g = G(0, 3, FWD);
		g.addEvidence(SCE(FWD, withQual(new byte[] { 4,4,4,4,4,4,4}, withSequence("TAAAGTC", Read(0, 10, "4M3S")))));
		g.addEvidence(SCE(FWD, withQual(new byte[] { 4,4,4,4,4,4,4}, withSequence("AAAGTCT", Read(0, 10, "3M4S")))));
		List<AssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		AssemblyEvidence bp = result.get(0);
		assertEquals(4, bp.getAssemblySoftClipLengthMax(null));
		assertEquals(7, bp.getAssemblySoftClipLengthTotal(null));
		assertEquals(2, bp.getAssemblySupportCountSoftClip(null));
	}
	public void should_set_assembly_mate_attributes() {
		DeBruijnReadGraph g = G(0, 3, FWD);
		g.addEvidence(NRRP(OEA(0, 1, "5M", true)));
		List<AssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		AssemblyEvidence bp = result.get(0);
		assertEquals(5, bp.getAssemblyReadPairLengthMax(null));
	}
	@Test
	public void should_assemble_adjacent_scs() {
		DeBruijnReadGraph g = G(0, 3, FWD);
		g.addEvidence(SCE(FWD, withSequence("AAAGTC", Read(0, 10, "3M3S")))); // 12
		g.addEvidence(SCE(FWD, withSequence("AAAAGTCTT", Read(0, 10, "4M4S")))); // 13
		g.addEvidence(SCE(FWD, withSequence("AAAAGTCTT", Read(0, 10, "4M4S")))); // 13
		List<AssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		assertEquals(1, result.size());
		AssemblyEvidence bp = result.get(0);
		assertEquals("AAAGTCTT", S(bp.getAssemblySequence()));
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
		List<AssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		// we should call to 12 since we have many reads supporting that assembly
		assertEquals(12, result.get(0).getBreakendSummary().start);
	}
	@Test
	public void base_quals_should_match_positive_strand() {
		DeBruijnReadGraph g = G(0, 3, BWD);
		g.addEvidence(SCE(BWD, withSequence(  "AAGTCTTT", new byte[]       { 1, 2, 3, 4, 5, 6, 7, 8}, Read(0, 10, "5S3M"))));
		g.addEvidence(SCE(BWD, withSequence( "TAAGTCTTT", new byte[]    { 1, 2, 3, 4, 5, 6, 7, 8, 9}, Read(0, 10, "6S3M"))));
		g.addEvidence(SCE(BWD, withSequence("GTAAGTCTTT", new byte[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, Read(0, 10, "7S3M"))));
		List<AssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		assertEquals("GTAAGTC", S(result.get(0).getBreakendSequence()));
		assertEquals("GTAAGTCTTT", S(result.get(0).getAssemblySequence()));
		assertArrayEquals(new byte[] { 1, 1, 1, 3, 6, 9, 12, }, result.get(0).getBreakendQuality());
		
		g = G(0, 3, FWD);
		g.addEvidence(SCE(FWD, withSequence("TTTCTGAA",  new byte[]      { 8, 7, 6, 5, 4, 3, 2, 1}, Read(0, 10, "3M5S"))));
		g.addEvidence(SCE(FWD, withSequence("TTTCTGAAT", new byte[]    {9, 8, 7, 6, 5, 4, 3, 2, 1}, Read(0, 10, "3M6S"))));
		g.addEvidence(SCE(FWD, withSequence("TTTCTGAATG",new byte[] {10,9, 8, 7, 6, 5, 4, 3, 2, 1}, Read(0, 10, "3M7S"))));
		result = Lists.newArrayList(g.assembleContigsBefore(10000));
		assertEquals("CTGAATG", S(result.get(0).getBreakendSequence()));
		assertEquals("TTTCTGAATG", S(result.get(0).getAssemblySequence()));
		assertArrayEquals(new byte[] { 12, 9, 6, 3, 1, 1, 1 }, result.get(0).getBreakendQuality());
	}
	@Test
	public void assembly_count_should_match_expected() {
		DeBruijnReadGraph g = G(0, 3, BWD);
		g.addEvidence(SCE(BWD, withSequence(  "AAGTCTTT", Read(0, 10, "5S3M"))));
		g.addEvidence(SCE(BWD, withSequence( "TAAGTCTTT", Read(0, 10, "6S3M"))));
		g.addEvidence(SCE(BWD, withSequence("GTAAGTCTTT", Read(0, 10, "7S3M"))));
		List<AssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		assertEquals(5 + 6 + 7, result.get(0).getAssemblyBaseCount(EvidenceSubset.ALL));
		assertEquals(3, result.get(0).getAssemblySupportCountSoftClip(EvidenceSubset.ALL));
		assertEquals(3, result.get(0).getAssemblyAnchorLength());
	}
	@Test
	public void should_assemble_sc_and_nrrp_together() {
		DeBruijnReadGraph g = G(0, 4, FWD);
		g.addEvidence(SCE(FWD, withSequence("TAAAGTCC", Read(0, 10, "4M4S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAAGTCCT", Read(0, 10, "4M5S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAAGTCCT", Read(0, 10, "4M5S"))));
		g.addEvidence(NRRP(withSequence("GTCCTAGAC", DP(0, 1, "8M", true, 1, 10, "8M", false))));
		List<AssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		assertEquals(1, result.size());
		AssemblyEvidence bp = result.get(0);
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
		assertEquals("TAAAAGTCCTAGAC", S(bp.getAssemblySequence()));
	}
	@Test
	public void should_assemble_nrrp() {
		DeBruijnReadGraph g = G(0, 3, FWD);
		// expect FR orientation so we need to reverse comp the unmapped seq
		g.addEvidence(NRRP(withSequence("AGAC", OEA(0, 10, "4M", true)))); // revcomp=GTCT
		g.addEvidence(NRRP(withSequence("CTAG", DP(0, 9, "4M", true, 1, 10, "4M", false))));
		List<AssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		assertEquals(1, result.size());
		AssemblyEvidence bp = result.get(0);
		assertEquals("GTCTAG", S(bp.getAssemblySequence()));
		assertEquals(FWD, bp.getBreakendSummary().direction);
		assertEquals(0, bp.getBreakendSummary().referenceIndex);
		assertEquals(13, bp.getBreakendSummary().start);
	}
	@Test
	public void should_assemble_best_contig() {
		DeBruijnReadGraph g = G(0, 3, FWD);
		g.addEvidence(SCE(FWD, withSequence("AAAGTCTA", Read(0, 10, "3M5S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAGTCTA", Read(0, 10, "3M5S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAGTCTG", Read(0, 10, "3M5S"))));
		List<AssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		assertEquals(1, result.size());
		AssemblyEvidence bp = result.get(0);
		// CTG becomes unanchored & we can't make a contig out of it
		assertEquals("AAAGTCTA", S(bp.getAssemblySequence()));
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
		List<AssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		assertEquals(1, result.size()); // Update: we now remove the entire tree when do the first assembly
		
		//assertEquals(2, result.size());
		//assertEquals("TAAAAGTCCTAGAC", result.get(1).getAssemblyConsensus());
		// AGAT is left over as a seed but since it does not anchor to the reference
		// As it's shorter than a read assembly, we presume that it is
		// sequencing noise and abandon it
		//assertEquals("AGAT", result.get(0).getAssemblyConsensus());
		//assertTrue(result.get(0).getFilters().contains(VcfFilter.ASSEMBLY_TOO_SHORT.name()));
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
		List<AssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		assertEquals(1, result.size());
		AssemblyEvidence bp = result.get(0);
		// CTG becomes unanchored & we can't make a contig out of it
		assertEquals(1, bp.getAssemblySupportCountSoftClip(null));
	}
	@Test
	public void should_exclude_sc_nrrp_without_ref_anchor() {
		DeBruijnReadGraph g = G(0, 3, FWD);
		// main assembly
		g.addEvidence(SCE(FWD, withSequence("AAAAGTCCTAG", Read(0, 10, "4M7S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAAGTCCTAG", Read(0, 10, "4M7S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAAGTCCTAG", Read(0, 10, "4M7S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAAGTCCTAG", Read(0, 10, "4M7S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAAGTCCTAG", Read(0, 10, "4M7S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAAGTCCTAG", Read(0, 10, "4M7S"))));
		// problematic assembly: we have an anchor sequence, but we've used it elsewhere
		// so we end up unanchored
		g.addEvidence(SCE(FWD, withSequence("AAAAGTCCTATG", Read(0, 10, "4M7S"))));
		g.addEvidence(NRRP(withSequence(             "ATGT", DP(0, 1, "8M", true, 1, 10, "8M", false))));
		List<AssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		assertEquals(1, result.size());
	}
	@Test
	public void should_filter_sc_branch() {
		DeBruijnReadGraph g = G(0, 3, FWD);
		// main assembly
		g.addEvidence(SCE(FWD, withSequence("AAAAGTCCTAG", Read(0, 10, "4M7S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAAGTCCTAG", Read(0, 10, "4M7S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAAGTCCTAG", Read(0, 10, "4M7S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAAGTCCTAG", Read(0, 10, "4M7S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAAGTCCTAG", Read(0, 10, "4M7S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAAGTCCTAG", Read(0, 10, "4M7S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAAGTCCTATG", Read(0, 10, "4M7S"))));
		List<AssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		assertEquals(1, result.size());
		// TATG should not be included as a result 
	}
	@Test
	public void removeBefore_should_remove_kmers() {
		DeBruijnReadGraph g = G(0, 3, FWD);
		g.addEvidence(SCE(FWD, withSequence("AAAGTC", Read(0, 9, "3M3S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAGTC", Read(0, 10, "3M3S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAGTC", Read(0, 10, "3M3S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAGTC", Read(0, 10, "3M3S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAGTC", Read(0, 11, "3M3S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAAGTCT", Read(0, 10, "4M4S"))));
		g.assembleContigsBefore(10000);
		g.removeBefore(10000);
		assertEquals(0, g.size());
	}
	@Test
	public void should_count_only_breakend_evidence() {
		DeBruijnReadGraph g = G(0, 3, FWD);
		g.addEvidence(SCE(FWD, SES(true),  withSequence("TATG", Read(0, 10, "3M1S"))));
		g.addEvidence(SCE(FWD, SES(false), withSequence("TATT", Read(0, 10, "3M1S"))));
		List<AssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(10000));
		assertEquals(2, result.size());
		assertEquals(1, result.get(0).getAssemblySupportCountSoftClip(EvidenceSubset.ALL));
		assertEquals(1, result.get(1).getAssemblySupportCountSoftClip(EvidenceSubset.ALL));
	}
	@Test
	public void removeBefore_should_clear_unanchored_subgraphs_SC_N_bases() {
		DeBruijnReadGraph g = G(0, 14, FWD);
		g.addEvidence(SCE(FWD, SES(true),  withSequence("CATTAATCGCAAGAGCGGGTTGTATTCGACNNNCAAGTCAGCTGAAGCACCATTACCCGATCANAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(0, 10, "15M85S"))));
		g.removeBefore(0);
		g.sanityCheckSubgraphs();
	}
	@Test
	public void subgraphs_must_be_anchored_RP_N_bases() {
		DeBruijnReadGraph g = G(0, 14, FWD);
		g.addEvidence(NRRP(withSequence("CATTAATCGCAAGAGCGGGTTGTATTCGACNNNCAAGTCAGCTGAAGCACCATTACCCGATCANAACATATCAGAAATGATTGACGTATCACAAGCCGGA", OEA(0, 1, "100M", true))));
		g.sanityCheckSubgraphs();
	}
}
