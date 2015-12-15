package au.edu.wehi.idsv.debruijn.subgraph;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.AssemblyEvidence;
import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.SAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.SoftClipEvidence;
import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.configuration.AssemblyConfiguration;

import com.google.common.collect.Lists;


public class DeBruijnReadGraphTest extends TestHelper {
	public DeBruijnReadGraph G(int referenceIndex, int k) {
		AssemblyEvidenceSource aes = AES();
		AssemblyConfiguration p = aes.getContext().getConfig().getAssembly();
		p.errorCorrection.maxBaseMismatchForCollapse = 0;
		p.subgraph.traversalMaximumBranchingFactor = 1;
		p.k = k;
		return new DeBruijnReadGraph(aes, referenceIndex, null);
	}
	@Test
	public void should_exclude_kmers_containing_ambiguous_base() {
		DeBruijnReadGraph g = G(0, 4);
		g.addEvidence(SCE(FWD, withSequence("TANAGTN", Read(0, 10, "4M3S"))));
		g.addEvidence(SCE(FWD, withSequence("AANNTCT", Read(0, 11, "3M4S"))));
		List<SAMRecordAssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(LCCB + 10000));
		assertEquals(0, result.size());
	}
	@Test
	public void should_assemble_sc() {
		DeBruijnReadGraph g = G(0, 3);
		g.addEvidence(SCE(FWD, withSequence("TAAAGTC", Read(0, 10, "4M3S"))));
		g.addEvidence(SCE(FWD, withSequence("AAAGTCT", Read(0, 11, "3M4S"))));
		List<SAMRecordAssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(LCCB + 10000));
		assertEquals(1, result.size());
		AssemblyEvidence bp = result.get(0);
		assertEquals("TAAAGTCT", S(bp.getAssemblySequence()));
	}
	@Test
	public void should_assemble_sc_b() {
		DeBruijnReadGraph g = G(0, 3);
		g.addEvidence(SCE(BWD, withSequence("TATG", Read(0, 10, "1S3M"))));
		g.addEvidence(SCE(BWD, withSequence("TTATG", Read(0, 10, "2S3M"))));
		List<SAMRecordAssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(LCCB + 10000));
		assertEquals(1, result.size());
		AssemblyEvidence bp = result.get(0);
		assertEquals("TTATG", S(bp.getAssemblySequence()));
	}
	@Test
	public void should_track_evidence() {
		List<SoftClipEvidence> evidence = Lists.newArrayList(
				SCE(FWD, withQual(new byte[] { 4,4,4,4,4,4,4}, withSequence("AAAGTCT", Read(0, 10, "3M4S")))),
				SCE(FWD, withQual(new byte[] { 4,4,4,4,4,4,4}, withSequence("TAAGTCT", Read(0, 10, "4M3S")))));
		DeBruijnReadGraph g = G(0, 3);
		for (SoftClipEvidence e : evidence) g.addEvidence(e);
		List<SAMRecordAssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(LCCB + 10000));
		AssemblyEvidence bp = result.get(0);
		for (SoftClipEvidence e : evidence) {
			assertTrue(bp.isPartOfAssembly(e));
		}
		assertFalse(bp.isPartOfAssembly(SCE(FWD, Read(0, 11, "4M4S"))));
	}
	@Test
	public void should_assemble_adjacent_scs() {
		DeBruijnReadGraph g = G(0, 3);
		g.addEvidence(SCE(FWD, withName("r1", withSequence("AAAGTC", Read(0, 10, "4M2S"))))); // 13
		g.addEvidence(SCE(FWD, withName("r2", withSequence("AAAAGTCTT", Read(0, 10, "5M4S"))))); // 14
		g.addEvidence(SCE(FWD, withName("r3", withSequence("AAAAGTCTT", Read(0, 10, "5M4S"))))); // 14
		List<SAMRecordAssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(LCCB + 10000));
		assertEquals(1, result.size());
		AssemblyEvidence bp = result.get(0);
		assertEquals("AAAGTCTT", S(bp.getAssemblySequence()));
		assertEquals(0, bp.getBreakendSummary().referenceIndex);
		assertEquals(14, bp.getBreakendSummary().start);
		assertEquals(14, bp.getBreakendSummary().end);
	}
	@Test
	public void should_anchor_to_best_reference_position() {
		DeBruijnReadGraph g = G(0, 3);
		g.addEvidence(SCE(FWD, withName("r1", withSequence("AAAGTC", Read(0, 9, "3M3S")))));
		g.addEvidence(SCE(FWD, withName("r2", withSequence("AAAGTC", Read(0, 10, "3M3S")))));
		g.addEvidence(SCE(FWD, withName("r3", withSequence("AAAGTC", Read(0, 10, "3M3S")))));
		g.addEvidence(SCE(FWD, withName("r4", withSequence("AAAGTC", Read(0, 10, "3M3S")))));
		g.addEvidence(SCE(FWD, withName("r5", withSequence("AAAGTC", Read(0, 11, "3M3S")))));
		g.addEvidence(SCE(FWD, withName("r6", withSequence("AAAAGTCT", Read(0, 10, "4M4S")))));
		List<SAMRecordAssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(LCCB + 10000));
		// we should call to 12 since we have many reads supporting that assembly
		assertEquals(12, result.get(0).getBreakendSummary().start);
	}
	@Test
	public void base_quals_should_match_positive_strand() {
		DeBruijnReadGraph g = G(0, 3);
		g.addEvidence(SCE(BWD, withName("r1", withSequence(  "AAGTCTTT", new byte[]       { 1, 2, 3, 4, 5, 6, 7, 8}, Read(0, 10, "5S3M")))));
		g.addEvidence(SCE(BWD, withName("r2", withSequence( "TAAGTCTTT", new byte[]    { 1, 2, 3, 4, 5, 6, 7, 8, 9}, Read(0, 10, "6S3M")))));
		g.addEvidence(SCE(BWD, withName("r3", withSequence("GTAAGTCTTT", new byte[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, Read(0, 10, "7S3M")))));
		List<SAMRecordAssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(LCCB + 10000));
		assertEquals("GTAAGTC", S(result.get(0).getBreakendSequence()));
		assertEquals("GTAAGTCTTT", S(result.get(0).getAssemblySequence()));
		//     1 2  3  4  5  6
		//   1 2 3  4  5  6  7
		// 1 2 3 4  5  6  7  8
		// ---------------
		// 1 3 6 9 12 15 18 21 kmer weights
		// 
		assertArrayEquals(new byte[] {
				1,
				(1+3)/2,
				(1+3+6)/3,
				(3+6+9)/3,
				(6+9+12)/3,
				(9+12+15)/3,
				(12+15+18)/3/*,
				(15+18+21)/3,
				(18+21)/2,
				21*/
				}, result.get(0).getBreakendQuality());
		
		g = G(0, 3);
		g.addEvidence(SCE(FWD, withName("r1", withSequence("TTTCTGAA",  new byte[] {8, 7, 6, 5, 4, 3, 2, 1},       Read(0, 10, "3M5S")))));
		g.addEvidence(SCE(FWD, withName("r2", withSequence("TTTCTGAAT", new byte[] {9, 8, 7, 6, 5, 4, 3, 2, 1},    Read(0, 10, "3M6S")))));
		g.addEvidence(SCE(FWD, withName("r3", withSequence("TTTCTGAATG",new byte[] {10,9, 8, 7, 6, 5, 4, 3, 2, 1}, Read(0, 10, "3M7S")))));
		result = Lists.newArrayList(g.assembleContigsBefore(LCCB + 10000));
		assertEquals("CTGAATG", S(result.get(0).getBreakendSequence()));
		assertEquals("TTTCTGAATG", S(result.get(0).getAssemblySequence()));
		// 6  5  4  3  2  1
		// 7  6  5  4  3  2  1
		// 8  7  6  5  4  3  2  1
		//-----------------------
		//21 18 15 12  9  6  3  1
		
		assertArrayEquals(new byte[] {
				(18+15+12)/3,
				(15+12+9)/3,
				(12+9+6)/3,
				(9+6+3)/3,
				(6+3+1)/3,
				(3+1)/2,
				1
				}, result.get(0).getBreakendQuality());
	}
	@Test
	public void assembly_count_should_match_expected() {
		DeBruijnReadGraph g = G(0, 3);
		List<DirectedEvidence> evidence = Lists.newArrayList();
		evidence.add(SCE(BWD, withName("r1", withSequence(  "AAGTCTTT", Read(0, 10, "5S3M")))));
		evidence.add(SCE(BWD, withName("r2", withSequence( "TAAGTCTTT", Read(0, 10, "6S3M")))));
		evidence.add(SCE(BWD, withName("r3", withSequence("GTAAGTCTTT", Read(0, 10, "7S3M")))));
		for (DirectedEvidence e : evidence) g.addEvidence(e);
		List<SAMRecordAssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(LCCB + 10000));
		result.get(0).hydrateEvidenceSet(evidence).annotateAssembly();
		assertEquals(3, result.get(0).getAssemblySupportCountSoftClip());
		assertEquals(3, result.get(0).getAssemblyAnchorLength());
	}
	@Test
	public void should_assemble_sc_and_nrrp_together() {
		DeBruijnReadGraph g = G(0, 4);
		g.addEvidence(SCE(FWD, withName("r1", withSequence("TAAAGTCC", Read(0, 10, "4M4S")))));
		g.addEvidence(SCE(FWD, withName("r2", withSequence("AAAAGTCCT", Read(0, 10, "4M5S")))));
		g.addEvidence(SCE(FWD, withName("r3", withSequence("AAAAGTCCT", Read(0, 10, "4M5S")))));
		g.addEvidence(NRRP(withName("r4", withSequence("GTCCTAGAC", DP(0, 1, "9M", true, 1, 10, "9M", false)))));
		List<SAMRecordAssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(LCCB + 10000));
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
		DeBruijnReadGraph g = G(0, 3);
		// expect FR orientation so we need to reverse comp the unmapped seq
		g.addEvidence(NRRP(withName("r1", withSequence("AGAC", OEA(0, 10, "4M", true))))); // revcomp=GTCT
		g.addEvidence(NRRP(withName("r2", withSequence("CTAG", DP(0, 9, "4M", true, 1, 10, "4M", false)))));
		List<SAMRecordAssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(LCCB + 10000));
		assertEquals(1, result.size());
		AssemblyEvidence bp = result.get(0);
		assertEquals("GTCTAG", S(bp.getAssemblySequence()));
		assertEquals(FWD, bp.getBreakendSummary().direction);
		assertEquals(0, bp.getBreakendSummary().referenceIndex);
		assertEquals(13, bp.getBreakendSummary().start);
	}
	@Test
	public void should_assemble_best_contig() {
		DeBruijnReadGraph g = G(0, 3);
		g.addEvidence(SCE(FWD, withName("r1", withSequence("AAAGTCTA", Read(0, 10, "3M5S")))));
		g.addEvidence(SCE(FWD, withName("r2", withSequence("AAAGTCTA", Read(0, 10, "3M5S")))));
		g.addEvidence(SCE(FWD, withName("r3", withSequence("AAAGTCTG", Read(0, 10, "3M5S")))));
		List<SAMRecordAssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(LCCB + 10000));
		assertEquals(1, result.size());
		AssemblyEvidence bp = result.get(0);
		// CTG becomes unanchored & we can't make a contig out of it
		assertEquals("AAAGTCTA", S(bp.getAssemblySequence()));
	}
	@Test
	public void should_assemble_when_out_of_scope() {
		DeBruijnReadGraph g = G(0, 3);
		g.addEvidence(SCE(FWD, withName("r1", withSequence("TAAAGTC", Read(0, 1, "4M3S")))));
		g.addEvidence(SCE(FWD, withName("r2", withSequence("AAAGTCT", Read(0, 2, "3M4S")))));
		// reference support to kmer starting at ref position 6
		assertEquals(0, Lists.newArrayList(g.assembleContigsBefore(LCCB + 6)).size());
		assertEquals(1, Lists.newArrayList(g.assembleContigsBefore(LCCB + 7)).size());
	}
	@Test
	public void removeBefore_should_remove_entire_subgraph_from_debruijn_graph() {
		DeBruijnReadGraph g = G(0, 3);
		g.addEvidence(SCE(FWD, withName("r1", withSequence("TAAAGTC", Read(0, 1, "4M3S")))));
		g.addEvidence(SCE(FWD, withName("r2", withSequence("AAAGTCT", Read(0, 2, "3M4S")))));
		assertEquals(1, Lists.newArrayList(g.assembleContigsBefore(LCCB + 7)).size());
		g.removeBefore(LCCB + 7);
		assertEquals(0, Lists.newArrayList(g.assembleContigsBefore(LCCB + 7)).size());
		assertEquals(0, Lists.newArrayList(g.assembleContigsBefore(Long.MAX_VALUE)).size());
	}
	@Test
	public void should_filter_short_contigs_that_do_not_include_reference_bases() {
		DeBruijnReadGraph g = G(0, 4);
		g.addEvidence(SCE(FWD, withName("r1", withSequence("TAAAGTCC", Read(0, 10, "4M4S")))));
		g.addEvidence(SCE(FWD, withName("r2", withSequence("AAAAGTCCT", Read(0, 10, "4M5S")))));
		g.addEvidence(SCE(FWD, withName("r3", withSequence("AAAAGTCCT", Read(0, 10, "4M5S")))));
		g.addEvidence(NRRP(withName("r4", withSequence("GTCCTAGAC", DP(0, 1, "9M", true, 1, 10, "9M", false)))));
		g.addEvidence(NRRP(withName("r5", withSequence("GTCCTAGAC", DP(0, 1, "9M", true, 1, 10, "9M", false)))));
		g.addEvidence(NRRP(withName("r6", withSequence("GTCCTAGAC", DP(0, 1, "9M", true, 1, 10, "9M", false)))));
		g.addEvidence(NRRP(withName("r7", withSequence("GTCCTAGAT", DP(0, 1, "9M", true, 1, 10, "9M", false)))));
		g.addEvidence(NRRP(withName("r8", withSequence("GTCCTAGAT", DP(0, 1, "9M", true, 1, 10, "9M", false)))));
		List<SAMRecordAssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(LCCB + 10000));
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
		DeBruijnReadGraph g = G(0, 4);
		g.addEvidence(SCE(FWD, withName("r1", withSequence("TAAAGTCC", Read(0, 10, "4M4S")))));
		g.addEvidence(SCE(FWD, withName("r2", withSequence("AAAAGTCCT", Read(0, 10, "4M5S")))));
		g.addEvidence(SCE(FWD, withName("r3", withSequence("AAAAGTCCT", Read(0, 10, "4M5S")))));
		g.addEvidence(NRRP(withName("r4", withSequence("GTCCTAGAC", DP(0, 1, "9M", true, 1, 10, "9M", false)))));
		g.addEvidence(NRRP(withName("r5", withSequence("GTCCTAGAC", DP(0, 1, "9M", true, 1, 10, "9M", false)))));
		g.addEvidence(NRRP(withName("r6", withSequence("GTCCTAGAC", DP(0, 1, "9M", true, 1, 10, "9M", false)))));
		g.addEvidence(NRRP(withName("r7", withSequence("GTCCTAGAT", DP(0, 1, "9M", true, 1, 10, "9M", false)))));
		g.addEvidence(NRRP(withName("r8", withSequence("GTCCTAGAT", DP(0, 1, "9M", true, 1, 10, "9M", false)))));
		String s = g.debugPrintPaths();
		assertNotNull(s);
	}
	@Test
	public void reference_kmer_support_should_not_be_considered_read_support() {
		DeBruijnReadGraph g = G(0, 3);
		List<DirectedEvidence> evidence = Lists.newArrayList();
		evidence.add(SCE(FWD, withName("r1", withSequence("AAAAAAAA", Read(0, 10, "3M5S")))));
		evidence.add(SCE(FWD, withName("r2", withSequence("AAAGTCTA", Read(0, 10, "3M5S")))));
		for (DirectedEvidence e : evidence) g.addEvidence(e);
		List<SAMRecordAssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(LCCB + 10000));
		assertEquals(1, result.size());
		SAMRecordAssemblyEvidence bp = result.get(0).hydrateEvidenceSet(evidence).annotateAssembly();
		// CTG becomes unanchored & we can't make a contig out of it
		assertEquals(1, bp.getAssemblySupportCountSoftClip());
	}
	@Test
	public void should_exclude_sc_nrrp_without_ref_anchor() {
		DeBruijnReadGraph g = G(0, 3);
		// main assembly
		g.addEvidence(SCE(FWD, withName("r1", withSequence("AAAAGTCCTAG", Read(0, 10, "4M7S")))));
		g.addEvidence(SCE(FWD, withName("r2", withSequence("AAAAGTCCTAG", Read(0, 10, "4M7S")))));
		g.addEvidence(SCE(FWD, withName("r3", withSequence("AAAAGTCCTAG", Read(0, 10, "4M7S")))));
		g.addEvidence(SCE(FWD, withName("r4", withSequence("AAAAGTCCTAG", Read(0, 10, "4M7S")))));
		g.addEvidence(SCE(FWD, withName("r5", withSequence("AAAAGTCCTAG", Read(0, 10, "4M7S")))));
		g.addEvidence(SCE(FWD, withName("r6", withSequence("AAAAGTCCTAG", Read(0, 10, "4M7S")))));
		// problematic assembly: we have an anchor sequence, but we've used it elsewhere
		// so we end up unanchored
		g.addEvidence(SCE(FWD, withName("r7", withSequence("AAAAGTCCTATG", Read(0, 10, "4M8S")))));
		g.addEvidence(NRRP(withName("r8", withSequence(             "ATGT", DP(0, 1, "4M", true, 1, 10, "4M", false)))));
		List<SAMRecordAssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(LCCB + 10000));
		assertEquals(1, result.size());
	}
	@Test
	public void should_filter_sc_branch() {
		DeBruijnReadGraph g = G(0, 3);
		// main assembly
		g.addEvidence(SCE(FWD, withName("z1", withSequence("AAAAGTCCTAG", Read(0, 10, "4M7S")))));
		g.addEvidence(SCE(FWD, withName("z2", withSequence("AAAAGTCCTAG", Read(0, 10, "4M7S")))));
		g.addEvidence(SCE(FWD, withName("z3", withSequence("AAAAGTCCTAG", Read(0, 10, "4M7S")))));
		g.addEvidence(SCE(FWD, withName("z4", withSequence("AAAAGTCCTAG", Read(0, 10, "4M7S")))));
		g.addEvidence(SCE(FWD, withName("z5", withSequence("AAAAGTCCTAG", Read(0, 10, "4M7S")))));
		g.addEvidence(SCE(FWD, withName("z6", withSequence("AAAAGTCCTAG", Read(0, 10, "4M7S")))));
		g.addEvidence(SCE(FWD, withName("z7", withSequence("AAAAGTCCTATG", Read(0, 10, "4M8S")))));
		List<SAMRecordAssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(LCCB + 10000));
		assertEquals(1, result.size());
		// TATG should not be included as a result 
	}
	@Test
	public void removeBefore_should_remove_kmers() {
		DeBruijnReadGraph g = G(0, 3);
		g.addEvidence(SCE(FWD, withName("z1", withSequence("AAAGTC", Read(0, 9, "3M3S")))));
		g.addEvidence(SCE(FWD, withName("z2", withSequence("AAAGTC", Read(0, 10, "3M3S")))));
		g.addEvidence(SCE(FWD, withName("z3", withSequence("AAAGTC", Read(0, 10, "3M3S")))));
		g.addEvidence(SCE(FWD, withName("z4", withSequence("AAAGTC", Read(0, 10, "3M3S")))));
		g.addEvidence(SCE(FWD, withName("z5", withSequence("AAAGTC", Read(0, 11, "3M3S")))));
		g.addEvidence(SCE(FWD, withName("z6", withSequence("AAAAGTCT", Read(0, 10, "4M4S")))));
		g.assembleContigsBefore(LCCB + 10000);
		g.removeBefore(LCCB + 10000);
		assertEquals(0, g.size());
	}
	@Test
	public void should_count_only_breakend_evidence() {
		DeBruijnReadGraph g = G(0, 3);
		List<DirectedEvidence> evidence = Lists.newArrayList();
		evidence.add(SCE(FWD, SES(true),  withName("r1", withSequence("TATG", Read(0, 10, "3M1S")))));
		evidence.add(SCE(FWD, SES(false), withName("r2", withSequence("TATT", Read(0, 10, "3M1S")))));
		for (DirectedEvidence e : evidence) g.addEvidence(e);
		List<SAMRecordAssemblyEvidence> result = Lists.newArrayList(g.assembleContigsBefore(LCCB + 10000));
		assertEquals(2, result.size());
		result.get(0).hydrateEvidenceSet(evidence).annotateAssembly();
		result.get(1).hydrateEvidenceSet(evidence).annotateAssembly();
		assertEquals(1, result.get(0).getAssemblySupportCountSoftClip());
		assertEquals(1, result.get(1).getAssemblySupportCountSoftClip());
	}
	@Test
	public void removeBefore_should_clear_unanchored_subgraphs_SC_N_bases() {
		DeBruijnReadGraph g = G(0, 14);
		g.addEvidence(SCE(FWD, SES(true),  withSequence("CATTAATCGCAAGAGCGGGTTGTATTCGACNNNCAAGTCAGCTGAAGCACCATTACCCGATCANAACATATCAGAAATGATTGACGTATCACAAGCCGGA", Read(0, 10, "15M85S"))));
		g.removeBefore(0);
		g.sanityCheckSubgraphs();
	}
	@Test
	public void subgraphs_must_be_anchored_RP_N_bases() {
		DeBruijnReadGraph g = G(0, 14);
		g.addEvidence(NRRP(withSequence("CATTAATCGCAAGAGCGGGTTGTATTCGACNNNCAAGTCAGCTGAAGCACCATTACCCGATCANAACATATCAGAAATGATTGACGTATCACAAGCCGGA", OEA(0, 1, "100M", true))));
		g.sanityCheckSubgraphs();
	}
}
