package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import org.junit.Test;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.ReadEvidenceAssemblerUtil;
import au.edu.wehi.idsv.RealignedBreakpoint;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.vcf.VcfAttributes;


public class RealignedBreakpointTest extends TestHelper {
	@Test(expected=IllegalArgumentException.class)
	public void should_throw_if_realigned_unmapped() {
		VariantContextDirectedBreakpoint dba = ReadEvidenceAssemblerUtil.create(getContext(), "test", 0, 1, BreakendDirection.Backward, B("AT"), B("ATA"), 1, 1, 0);
		dba = new VariantContextDirectedBreakpoint(getContext(), new VariantContextBuilder(dba).make());
		new RealignedBreakpoint(dba.getBreakendSummary(), Unmapped(2));
	}
	@Test
	public void should_set_realign_evidence() {
		RealignedBreakpoint rbp = new RealignedBreakpoint(new BreakendSummary(0, FWD, 1, 1, null), withMapq(10, Read(0, 1, "5M"))[0]);
		assertEquals(5, rbp.getBreakpointSummary().evidence.get(VcfAttributes.REALIGN_MAX_LENGTH));
		assertEquals(5, rbp.getBreakpointSummary().evidence.get(VcfAttributes.REALIGN_TOTAL_LENGTH));
		//assertEquals(10, rbp.getBreakpointSummary().evidence.get(EvidenceAttributes.REALIGN_MAX_MAPQ));
		assertEquals(10, rbp.getBreakpointSummary().evidence.get(VcfAttributes.REALIGN_TOTAL_MAPQ));
	}
	public RealignedBreakpoint test_seq(String originalBreakpointSequence, String cigar, BreakendDirection direction, boolean alignNegativeStrand, String expectedUntemplatedSequence) {
		SAMRecord r = Read(0, 1, cigar);
		r.setReadBases(alignNegativeStrand ? B(SequenceUtil.reverseComplement(originalBreakpointSequence)): B(originalBreakpointSequence));
		r.setReadNegativeStrandFlag(alignNegativeStrand);
		RealignedBreakpoint rbp = new RealignedBreakpoint(new BreakendSummary(0, direction, 1, 1, null), r);
		assertEquals(expectedUntemplatedSequence, rbp.getInsertedSequence());
		return rbp;
	}
	@Test
	public void inserted_sequence_should_be_relative_to_original_realigned_sequence() {
		test_seq("ATTCNNAGC", "4S3M3S", FWD, false, "ATTC");
		test_seq("ATTCNNAGC", "4S3M3S", FWD, true, "ATT");
		test_seq("ATTCNNAGC", "4S3M3S", BWD, false, "AGC");
		test_seq("ATTCNNAGC", "4S3M3S", BWD, true, "NAGC");
	}
	@Test
	public void breakpoint_window_size_should_correspond_to_microhomology_length() {
		SAMRecord r = Read(0, 1, "10M5S");
		r.setReadBases(B("TTTTTAAAAAAAAAT"));
		SAMRecord realign = Read(0, 100, "5M");
		r.setReadBases(B("AAAAT"));
		SoftClipEvidence sce = new SoftClipEvidence(getContext(), FWD, r);
		RealignedBreakpoint rbp = new RealignedBreakpoint(sce.getBreakendSummary(), realign);
		// breakpoint could be anywhere in the poly A microhomology
		assertEquals(5, rbp.getBreakpointSummary().start);
		assertEquals(14, rbp.getBreakpointSummary().end);
		assertEquals(95, rbp.getBreakpointSummary().start2);
		assertEquals(104, rbp.getBreakpointSummary().end2);
	}
}
