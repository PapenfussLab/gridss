package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;

import org.junit.Test;

import au.edu.wehi.idsv.sam.SamTags;


public class SmallIndelSAMRecordAssemblyEvidenceTest extends TestHelper {
	private SmallIndelSAMRecordAssemblyEvidence create(int position, String cigar, String bases, BreakendDirection dir) {
		SAMRecord r = Read(0, position, cigar);
		r.setAttribute(SamTags.ASSEMBLY_DIRECTION, dir.toChar());
		r.setReadBases(B(bases));
		r.setBaseQualityString(bases);
		r.setMappingQuality(35);
		SmallIndelSAMRecordAssemblyEvidence e = new SmallIndelSAMRecordAssemblyEvidence(AES(), r);
		return e;
	}
	private void check_matches(int position, String cigar, String bases, BreakendDirection dir,
			int anchorPos, String anchorCigar, String anchorSeq,
			int realignPos, String realignCigar, String realignSeq) {
		SmallIndelSAMRecordAssemblyEvidence e = create(position, cigar, bases, dir);
		
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
		check_matches(2, "2M3D4M", "AANNNN", FWD, 2, "2M4S", "AANNNN", 7, "4M", "NNNN");
		check_matches(2, "2M3D4M", "AANNNN", BWD, 7, "2S4M", "AANNNN", 2, "2M", "AA");
		
		// 123---4567890
		//  MMiiiMMMM
		check_matches(2, "2M3I4M", "AATTTNNNN", FWD, 2, "2M7S", "AATTTNNNN", 4, "3S4M", "TTTNNNN");
		check_matches(2, "2M3I4M", "AATTTNNNN", BWD, 4, "5S4M", "AATTTNNNN", 2, "2M3S", "AATTT");
	}
	@Test
	public void should_split_read_on_largest_indel() {
		// 1234567----8901234567890
		//  MMdddMiiiiMMMMM
		
		check_matches(2, "2M3D1M4I5M", "AAATTTTNNNNN", FWD, 2, "2M3D1M9S", "AAATTTTNNNNN", 8, "4S5M", "TTTTNNNNN");
		check_matches(2, "2M3D1M4I5M", "AAATTTTNNNNN", BWD, 8, "7S5M", "AAATTTTNNNNN", 2, "2M3D1M4S", "AAATTTT");
	}
	@Test
	public void split_read_mapq_should_be_source_mapq() {
		SmallIndelSAMRecordAssemblyEvidence e = create(0, "1M1D1M", "NN", FWD);
		assertEquals(35, e.getSAMRecord().getMappingQuality());
		assertEquals(35, e.getRemoteSAMRecord().getMappingQuality());
	}
}
