package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;

import org.junit.Test;

import au.edu.wehi.idsv.sam.SamTags;


public class SmallIndelSAMRecordAssemblyEvidenceTest extends TestHelper {
	public static SmallIndelSAMRecordAssemblyEvidence create(int position, String cigar, String bases) {
		SAMRecord r = Read(0, position, cigar);
		r.setReadBases(B(bases));
		r.setBaseQualityString(bases);
		MockDirectedEvidence evidence = new MockDirectedEvidence(new BreakendSummary(0, FWD, 1, 2));
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCEID_UNCATEGORISED, evidence.getEvidenceID());
		SmallIndelSAMRecordAssemblyEvidence e = new SmallIndelSAMRecordAssemblyEvidence(AES(), r);
		e.hydrateEvidenceSet(evidence);
		e.annotateAssembly();
		return e;
	}
	private void check_matches(int position, String cigar, String bases,
			int anchorPos, String anchorCigar, String anchorSeq,
			int realignPos, String realignCigar, String realignSeq) {
		SmallIndelSAMRecordAssemblyEvidence e = create(position, cigar, bases);
		
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
	public void should_split_read_on_largest_indel() {
		// 1234567----8901234567890
		//  MMdddMiiiiMMMMM
		
		check_matches(2, "2M3D1M4I5M", "AAATTTTNNNNN", 2, "2M3D1M9S", "AAATTTTNNNNN", 8, "4S5M", "TTTTNNNNN");
	}
	@Test
	public void split_read_mapq_should_be_source_mapq() {
		SmallIndelSAMRecordAssemblyEvidence e = create(0, "1M1D1M", "NN");
		assertEquals(e.getBackingRecord().getMappingQuality(), e.getSAMRecord().getMappingQuality());
		assertEquals(e.getBackingRecord().getMappingQuality(), e.getRemoteSAMRecord().getMappingQuality());
	}
	@Test
	public void remote_mapq_should_match_local() {
		SmallIndelSAMRecordAssemblyEvidence e = create(0, "1M1D1M", "NN");
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
		SmallIndelSAMRecordAssemblyEvidence e = create(10, "1M2P2N2P1M", "NN");
		assertEquals(new BreakpointSummary(0, FWD, 10, 10, 0, BWD, 9, 9), e.getBreakendSummary());
	}
	@Test
	public void remote_breakend_SAMRecord_should_have_start_coordinate_of_high_breakend() {
		assertEquals(13, create(10, "1M2D1M", "NN").getRemoteSAMRecord().getAlignmentStart());
		assertEquals(9, create(10, "1M2P2N2P1M", "NN").getRemoteSAMRecord().getAlignmentStart());
	}
	@Test
	public void should_not_have_assembly_direction() {
		SmallIndelSAMRecordAssemblyEvidence e = create(0, "1M1D1M", "NN");
		assertNull(SAMRecordAssemblyEvidence.getBreakendDirection(e.getBackingRecord()));
	}
	@Test
	public void remote_should_not_be_considered_remote_evidence_since_assembly_spans_entire_breakpoint() {
		SmallIndelSAMRecordAssemblyEvidence e = create(0, "1M1D1M", "NN");
		assertFalse(e.asRemote() instanceof RemoteEvidence);
	}
	@Test
	public void getEvidenceID_should_prefix_with_local_breakend_direction() {
		SmallIndelSAMRecordAssemblyEvidence e = create(0, "1M1D1M", "NN");
		assertEquals('f', e.getEvidenceID().charAt(0));
		assertEquals('b', e.asRemote().getEvidenceID().charAt(0));
	}
	@Test
	public void isSpanningAssembly() {
		SmallIndelSAMRecordAssemblyEvidence e = create(0, "1M1D1M", "NN");
		assertTrue(e.isSpanningAssembly());
		assertTrue(e.asRemote().isSpanningAssembly());
		assertTrue(AssemblyFactory.hydrate(AES(), e.getBackingRecord()).isSpanningAssembly());
		assertTrue(AssemblyFactory.hydrate(AES(), e.getSAMRecord()).isSpanningAssembly());
		assertTrue(AssemblyFactory.incorporateRealignment(getContext(), AssemblyFactory.hydrate(AES(), e.getSAMRecord()), e.getRemoteSAMRecord()).isSpanningAssembly());
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
		SmallIndelSAMRecordAssemblyEvidence e = (SmallIndelSAMRecordAssemblyEvidence)AssemblyFactory.hydrate(AES(), r);
		assertEquals("136M224S", e.getSAMRecord().getCigarString());
		assertEquals(10603, e.getSAMRecord().getAlignmentStart());
		assertEquals("88S136M", e.getRemoteSAMRecord().getCigarString());
	}
}
