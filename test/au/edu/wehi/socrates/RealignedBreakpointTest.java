package au.edu.wehi.socrates;

import static org.junit.Assert.assertEquals;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import org.junit.Test;

import au.edu.wehi.socrates.vcf.VcfAttributes;


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
	public void test_seq(String originalBreakpointSequence, String cigar, BreakendDirection direction, boolean alignNegativeStrand, String expectedUntemplatedSequence) {
		SAMRecord r = Read(0, 1, cigar);
		r.setReadBases(alignNegativeStrand ? B(SequenceUtil.reverseComplement(originalBreakpointSequence)): B(originalBreakpointSequence));
		r.setReadNegativeStrandFlag(alignNegativeStrand);
		RealignedBreakpoint rbp = new RealignedBreakpoint(new BreakendSummary(0, direction, 1, 1, null), r);
		assertEquals(expectedUntemplatedSequence, rbp.getInsertedSequence());
	}
	@Test
	public void inserted_sequence_should_be_relative_to_original_realigned_sequence() {
		test_seq("ATTCNNAGC", "4S3M3S", FWD, false, "ATTC");
		test_seq("ATTCNNAGC", "4S3M3S", FWD, true, "ATT");
		test_seq("ATTCNNAGC", "4S3M3S", BWD, false, "AGC");
		test_seq("ATTCNNAGC", "4S3M3S", BWD, true, "NAGC");
	}
}
