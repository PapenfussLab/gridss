package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;

import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.sam.SamTags;

import com.google.common.collect.ImmutableList;


public class SpanningSAMRecordAssemblyEvidenceTest extends TestHelper {
	public static SpanningSAMRecordAssemblyEvidence create(int position, String cigar, String bases) {
		SAMRecord r = Read(0, position, cigar);
		assertEquals(bases.length(), r.getCigar().getReadLength());
		r.setReadBases(B(bases));
		r.setBaseQualityString(bases);
		MockDirectedEvidence evidence = new MockDirectedEvidence(new BreakendSummary(0, FWD, 1, 2));
		r.setAttribute(SamTags.EVIDENCEID, evidence.getEvidenceID());
		SAMRecordAssemblyEvidence e = new SAMRecordAssemblyEvidence(AES(), r, ImmutableList.of());
		e.hydrateEvidenceSet(evidence);
		e.annotateAssembly();
		e.getBackingRecord().setMappingQuality(50);
		e.getSAMRecord().setMappingQuality(50);
		e.getRemoteSAMRecord().setMappingQuality(50);
		return e.getSpannedIndels().get(0);
	}
	private void check_matches(int position, String cigar, String bases,
			int anchorPos, String anchorCigar, String anchorSeq,
			int realignPos, String realignCigar, String realignSeq) {
		SpanningSAMRecordAssemblyEvidence e = create(position, cigar, bases);
		
		assertEquals(anchorPos, e.getSAMRecord().getAlignmentStart());
		assertEquals(anchorCigar, e.getSAMRecord().getCigarString());
		assertEquals(anchorSeq, S(e.getSAMRecord().getReadBases()));
		assertEquals(bases, SAMUtils.phredToFastq(e.getSAMRecord().getBaseQualities()));
		
		assertEquals(realignPos, e.getRemoteSAMRecord().getAlignmentStart());
		assertEquals(realignCigar, e.getRemoteSAMRecord().getCigarString());
		assertEquals(realignSeq, S(e.getRemoteSAMRecord().getReadBases()));
		assertEquals(e.getRemoteSAMRecord().getReadBases().length, e.getRemoteSAMRecord().getBaseQualities().length);
	}
	@Test
	public void should_act_as_split_read_mapping() {
		// 1234567890
		//  MMdddMMMM
		check_matches(2, "2M3D4M", "AANNNN", 2, "2M4S", "AANNNN", 7, "4M", "NNNN");

		
		// 123---4567890
		//  MMiiiMMMM
		check_matches(2, "2M3I4M", "AATTTNNNN", 2, "2M7S", "AATTTNNNN", 4, "3S4M", "TTTNNNN");
	}
	@Test
	public void split_read_mapq_should_be_source_mapq() {
		SpanningSAMRecordAssemblyEvidence e = create(0, "1M1D1M", "NN");
		assertEquals(e.getBackingRecord().getMappingQuality(), e.getSAMRecord().getMappingQuality());
		assertEquals(e.getBackingRecord().getMappingQuality(), e.getRemoteSAMRecord().getMappingQuality());
	}
	@Test
	public void remote_mapq_should_match_local() {
		SpanningSAMRecordAssemblyEvidence e = create(0, "1M1D1M", "NN");
		assertEquals(e.getRemoteMapq(), e.getLocalMapq());
	}
	@Test
	public void should_consider_adjacent_indel_operands_as_single_event() {
		// 123 4567890
		//  MMidddMMMM
		check_matches(2, "2M1I3D4M", "AATNNNN", 2, "2M5S", "AATNNNN", 7, "1S4M", "TNNNN");
	}
	/**
	 * placeholder value required as SAM unable to represent such alignments
	 */
	@Test
	public void should_treat_xPxNxP_as_placeholder_for_alignment_to_earlier_position_x_bp_before_expected_position() {
		SpanningSAMRecordAssemblyEvidence e = create(10, "1M2P2N2P1M", "NN");
		assertEquals(new BreakpointSummary(0, FWD, 10, 10, 0, BWD, 9, 9), e.getBreakendSummary());
	}
	@Test
	public void remote_breakend_SAMRecord_should_have_start_coordinate_of_high_breakend() {
		assertEquals(13, create(10, "1M2D1M", "NN").getRemoteSAMRecord().getAlignmentStart());
		assertEquals(9, create(10, "1M2P2N2P1M", "NN").getRemoteSAMRecord().getAlignmentStart());
	}
	@Test
	public void should_not_have_assembly_direction() {
		SpanningSAMRecordAssemblyEvidence e = create(0, "1M1D1M", "NN");
		assertEquals(FWD, e.getBreakendSummary().direction);
	}
	@Test
	public void remote_should_not_be_considered_remote_evidence_since_assembly_spans_entire_breakpoint() {
		SpanningSAMRecordAssemblyEvidence e = create(0, "1M1D1M", "NN");
		assertFalse(e.asRemote() instanceof RemoteEvidence);
	}
	@Test
	public void getEvidenceID_should_suffix_with_local_breakend_direction_and_indel_offset() {
		SpanningSAMRecordAssemblyEvidence e = create(0, "1M1D1M", "NN");
		String evidenceID = e.getParentAssembly().getEvidenceID();
		assertEquals(evidenceID + "_f0", e.getEvidenceID());
		assertEquals(evidenceID + "_b0", e.asRemote().getEvidenceID());
	}
	@Test
	public void qual_should_be_split_over_both_directions() {
		SoftClipEvidence sce = SCE(FWD, Read(0, 1, "1M5S"));
		SAMRecordAssemblyEvidence ass = AssemblyFactory.createAnchoredBreakpoint(getContext(), AES(), ImmutableList.of(sce.getEvidenceID()),
				0, 1, 1,
				0, 3, 1,
				B("AAA"), B("AAA"));
		ass.hydrateEvidenceSet(sce);
		ass.annotateAssembly();
		assertEquals(sce.getBreakendQual() / 2, ass.getSpannedIndels().get(0).getBreakendQual(), 0.0001);
		assertEquals(sce.getBreakendQual() / 2, ass.getSpannedIndels().get(0).getBreakpointQual(), 0.0001);
	}
	@Test
	public void read_pair_conversion_soft_clip_in_correct_direction() {
		// asm49752266     0       MT      10603   0       136M88I8952P8952N8952P136M      *       0       0       CTCTCATAACCCTCAACACCCACTCCCTCTTAGCCAATATTGTGCCTATTGCCATACTAGTCTTTGCCGCCTGCGAAGCAGCGGTGGGCCTAGCCCTACTAGTCTCAATCTCCAACACATATGGCCTAGACTACGGCTTAGTTAAACTTTCGTTTATTGCTAAAGGTTAATCACTGCTGTTTCCCGTGGGGGTGTGGCTAGGCTAAGCGTTTTGAGCTGCATTGCCAAGGGAAAGATGAAAAATTATAACCAAGCATAATATAGCAAGGACTAACCCCTATACCTTCTGCATAATGAATTAACTAGAAATAACTTTGCAAGGAGAGCCAAAGCTAAGACCCCCGAAACCAGACGAGCTAC        _______________________________________________________________________________________________________________________________________@<<<<<<<<<>ACEGIKMOQSUWY[]_______________________________________________________________________________________________________________________________________________________________________________________________________        bc:B:i,147      es:Z:RbST-E00106:108:H03M0ALXX:1:1118:24253:3735/2 fST-E00106:108:H03M0ALXX:1:1205:2583:30439/1 fST-E00106:108:H03M0ALXX:1:2207:15443:7989/2
		SAMRecord r = new SAMRecord(getContext().getBasicSamHeader());
		r.setReadName("asm49752266");
		r.setAlignmentStart(10603);
		r.setReferenceIndex(0);
		r.setCigarString("136M88I8952P8952N8952P136M");
		r.setReadBases(B("CTCTCATAACCCTCAACACCCACTCCCTCTTAGCCAATATTGTGCCTATTGCCATACTAGTCTTTGCCGCCTGCGAAGCAGCGGTGGGCCTAGCCCTACTAGTCTCAATCTCCAACACATATGGCCTAGACTACGGCTTAGTTAAACTTTCGTTTATTGCTAAAGGTTAATCACTGCTGTTTCCCGTGGGGGTGTGGCTAGGCTAAGCGTTTTGAGCTGCATTGCCAAGGGAAAGATGAAAAATTATAACCAAGCATAATATAGCAAGGACTAACCCCTATACCTTCTGCATAATGAATTAACTAGAAATAACTTTGCAAGGAGAGCCAAAGCTAAGACCCCCGAAACCAGACGAGCTAC"));
		r.setBaseQualities(B("_______________________________________________________________________________________________________________________________________@<<<<<<<<<>ACEGIKMOQSUWY[]_______________________________________________________________________________________________________________________________________________________________________________________________________"));
		r.setAttribute("bc", 147);
		r.setAttribute("es", "RbST-E00106:108:H03M0ALXX:1:1118:24253:3735/2 fST-E00106:108:H03M0ALXX:1:1205:2583:30439/1 fST-E00106:108:H03M0ALXX:1:2207:15443:7989/2");
		r.setMappingQuality(50);
		List<SpanningSAMRecordAssemblyEvidence> indels = AssemblyFactory.hydrate(AES(), r).getSpannedIndels();
		assertEquals(2, indels.size());
		SpanningSAMRecordAssemblyEvidence e = indels.get(0);
		assertEquals("136M224S", e.getSAMRecord().getCigarString());
		assertEquals(10603, e.getSAMRecord().getAlignmentStart());
		assertEquals("88S136M", e.getRemoteSAMRecord().getCigarString());
	}
	@Test
	public void indels_should_have_different_names() {
		SpanningSAMRecordAssemblyEvidence e = create(0, "1M1D1M1I1M", "NNNN");
		assertEquals(4, e.getParentAssembly().getSpannedIndels().size());
		assertEquals(4, e.getParentAssembly().getSpannedIndels().stream().map(x -> x.getEvidenceID()).distinct().count());
	}
}
