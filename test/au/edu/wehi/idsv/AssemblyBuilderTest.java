package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import au.edu.wehi.idsv.vcf.VcfAttributes;
import au.edu.wehi.idsv.vcf.VcfSvConstants;


public class AssemblyBuilderTest extends TestHelper {
	@Test
	public void should_filter_fully_reference_assemblies() {
		AssemblyBuilder ab = new AssemblyBuilder(getContext(), AES())
			.assemblyBases(B("AA"))
			.anchorLength(2)
			.referenceAnchor(0, 1)
			.direction(BWD);
		assertTrue(ab.makeVariant().isFiltered());
	}
	@Test
	public void should_filter_single_read_assemblies() {
		AssemblyBuilder ab = new AssemblyBuilder(getContext(), AES())
			.assemblyBases(B("AA"))
			.anchorLength(1)
			.referenceAnchor(0, 1)
			.direction(BWD)
			.assembledReadCount(1);
		assertTrue(ab.makeVariant().isFiltered());
	}
	@Test
	public void should_filter_mate_anchored_assembly_shorter_than_read_length() {
		AssemblyBuilder ab = new AssemblyBuilder(getContext(), AES())
			.assemblyBases(B("AA"))
			.mateAnchor(0, 1)
			.direction(BWD)
			.longestSupportingRead(3);
		assertTrue(ab.makeVariant().isFiltered());
	}
	//@Test // don't filter on assembly size
	public void should_filter_breakend_assembly_no_longer_than_soft_clip_size() {
		AssemblyBuilder ab = new AssemblyBuilder(getContext(), AES())
			.assemblyBases(B("AA"))
			.anchorLength(1)
			.referenceAnchor(0, 1)
			.direction(BWD)
			.maximumSoftClipLength(1);
		assertTrue(ab.makeVariant().isFiltered());
	}
	@Test
	public void read_length_filter_should_not_apply_to_anchored_breakend_assembly() {
		AssemblyBuilder ab = new AssemblyBuilder(getContext(), AES())
			.assemblyBases(B("AA"))
			.anchorLength(1)
			.referenceAnchor(0, 1)
			.direction(BWD)
			.longestSupportingRead(2);
		assertFalse(ab.makeVariant().isFiltered());
	}
	@Test
	public void soft_clip_size_filter_should_not_apply_to_unanchored_assembly() {
		AssemblyBuilder ab = new AssemblyBuilder(getContext(), AES())
			.assemblyBases(B("AA"))
			.mateAnchor(0, 1)
			.direction(BWD)
			.maximumSoftClipLength(1);
		assertFalse(ab.makeVariant().isFiltered());
	}
	@Test
	public void should_assembly_base_qualities() {
		AssemblyBuilder ab = new AssemblyBuilder(getContext(), AES())
			.assemblyBases(B("AA"))
			.referenceAnchor(0, 1)
			.direction(BWD);
		assertFalse(ab.makeVariant().isFiltered());
	}
	@Test
	public void should_set_assembly_evidence() {
		VariantContextDirectedBreakpoint dba = new AssemblyBuilder(getContext(), AES())
			.assemblyBases(B("AAAAAA"))
			.referenceAnchor(0, 1)
			.direction(BWD)
			.assembledReadCount(4)
			.assembledBaseCount(5)
			.makeVariant();
		assertEquals(4, dba.getBreakendSummary().evidence.get(VcfAttributes.ASSEMBLY_READS));
		assertEquals(5, dba.getBreakendSummary().evidence.get(VcfAttributes.ASSEMBLY_BASES));
		assertEquals(6, dba.getBreakendSummary().evidence.get(VcfAttributes.ASSEMBLY_LENGTH));		
	}
	@Test
	public void should_set_assembly_consensus() {
		VariantContextDirectedBreakpoint dba = new AssemblyBuilder(getContext(), AES())
			.assemblyBases(B("GTAC"))
			.referenceAnchor(0, 1)
			.direction(BWD)
			.assembledReadCount(4)
			.assembledBaseCount(5)
			.makeVariant();
		assertEquals("GTAC", dba.getAssemblyConsensus());
	}
	@Test
	public void should_set_assembler_name() {
		VariantContextDirectedBreakpoint dba = new AssemblyBuilder(getContext(), AES())
			.assemblyBases(B("AAAAAA"))
			.referenceAnchor(0, 1)
			.direction(BWD)
			.assemblerName("testAssembler")
			.makeVariant();
		assertEquals("testAssembler", dba.getAssemblerProgram());
	}
	public void should_set_max_soft_clip_length() {
		VariantContextDirectedBreakpoint dba = new AssemblyBuilder(getContext(), AES())
			.assemblyBases(B("AAAAAA"))
			.referenceAnchor(0, 1)
			.direction(BWD)
			.maximumSoftClipLength(12)
			.makeVariant();
		assertEquals(12, dba.getAssemblyMaximumSoftClipLength());
	}
	public void should_set_longest_read_length() {
		VariantContextDirectedBreakpoint dba = new AssemblyBuilder(getContext(), AES())
			.assemblyBases(B("AAAAAA"))
			.mateAnchor(0, 1)
			.direction(BWD)
			.longestSupportingRead(11)
			.makeVariant();
		assertEquals(11, dba.getAssemblyLongestSupportingRead());
	}
	@Test
	public void assembly_quality_should_default_to_average_breakend_base_quality() {
		VariantContextDirectedBreakpoint dba = new AssemblyBuilder(getContext(), AES())
			.assemblyBases(B("AAAAAAAA"))
			.referenceAnchor(0, 10)
			.direction(FWD)
			.anchorLength(3)
			.assemblyBaseQuality(new byte[] { 0,1,2,3,4,5,6,7 })
			.makeVariant();
		assertEquals(5d, dba.getAssemblyQuality(), 0);
	}
	@Test
	public void should_set_breakend_sequence() {
		assertEquals("GT", new AssemblyBuilder(getContext(), AES())
			.assemblyBases(B("GTAC"))
			.referenceAnchor(0, 1)
			.direction(BWD)
			.anchorLength(2)
			.makeVariant().getBreakpointSequenceString());
		assertEquals("AC", new AssemblyBuilder(getContext(), AES())
			.assemblyBases(B("GTAC"))
			.referenceAnchor(0, 10)
			.direction(FWD)
			.anchorLength(2)
			.makeVariant().getBreakpointSequenceString());
	}
	@Test
	public void should_breakend_to_reference_anchor() {
		VariantContextDirectedBreakpoint dba;
		BreakendSummary bs;
		dba = new AssemblyBuilder(getContext(), AES())
			.assemblyBases(B("AAAAAA"))
			.referenceAnchor(1, 10)
			.direction(BWD)
			.anchorLength(2)
			.makeVariant();
		bs = dba.getBreakendSummary();
		assertEquals(BWD, bs.direction);
		assertEquals(1, bs.referenceIndex);
		assertEquals(10, bs.start);
		assertEquals(10, bs.end);
		
		dba = new AssemblyBuilder(getContext(), AES())
			.assemblyBases(B("AAAAAA"))
			.referenceAnchor(1, 10)
			.direction(FWD)
			.anchorLength(2)
			.makeVariant();
		bs = dba.getBreakendSummary();
		assertEquals(FWD, bs.direction);
		assertEquals(1, bs.referenceIndex);
		assertEquals(10, bs.start);
		assertEquals(10, bs.end);
	}
	@Test
	public void should_set_breakend_to_mate_anchor_interval() {
		VariantContextDirectedBreakpoint dba;
		BreakendSummary bs;
		dba = new AssemblyBuilder(getContext(), AES())
			.assemblyBases(B("AAAAAA"))
			.mateAnchor(1, 1000)
			.direction(BWD)
			.makeVariant();
		bs = dba.getBreakendSummary();
		assertEquals(BWD, bs.direction);
		assertEquals(1, bs.referenceIndex);
		assertEquals(1000 - AES().getMetrics().getMaxFragmentSize(), bs.start);
		assertEquals(1000, bs.end);
		
		dba = new AssemblyBuilder(getContext(), AES())
			.assemblyBases(B("AAAAAA"))
			.mateAnchor(1, 1000)
			.direction(FWD)
			.makeVariant();
		bs = dba.getBreakendSummary();
		assertEquals(FWD, bs.direction);
		assertEquals(1, bs.referenceIndex);
		assertEquals(1000, bs.start);
		assertEquals(1000 + AES().getMetrics().getMaxFragmentSize(), bs.end);
	}
	@Test
	public void should_restrict_mate_anchor_interval_based_on_anchor_positions() {
		// max fragment size = 300
		VariantContextDirectedBreakpoint dba;
		BreakendSummary bs;
		dba = new AssemblyBuilder(getContext(), AES())
			.assemblyBases(B("AAAAAA"))
			.mateAnchor(1, 100, 5, 10, 20)
			.direction(FWD)
			.makeVariant();
		bs = dba.getBreakendSummary();
		fail(); // need to design this properly 
		//assertEquals(110, bs.start);
		//assertEquals(120, bs.end);
	}
	@Test
	public void mate_anchor_should_set_imprecise_header() {
		// max fragment size = 300
		VariantContextDirectedBreakpoint dba;
		BreakendSummary bs;
		dba = new AssemblyBuilder(getContext(), AES())
			.assemblyBases(B("AAAAAA"))
			.mateAnchor(1, 100)
			.direction(FWD)
			.makeVariant();
		bs = dba.getBreakendSummary();
		assertTrue(dba.hasAttribute(VcfSvConstants.IMPRECISE_KEY));
	}
}
