package au.edu.wehi.idsv;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import org.junit.Test;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.ReadEvidenceAssemblerUtil;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.vcf.VcfAttributes;


public class ReadEvidenceAssemblerUtilTest extends TestHelper {
	@Test
	public void anchor_should_use_reference_base_not_assembly_base() {
		String alt = ReadEvidenceAssemblerUtil.create(getContext(), "test", 0, 1, BreakendDirection.Forward, B("A"), B("CC"), 1, 1, 0).getAlternateAllele(0).getDisplayString();
		assertEquals('A', alt.charAt(0));
	}
	@Test
	public void should_assembly_base_qualities() {
		VariantContextDirectedBreakpoint dba = ReadEvidenceAssemblerUtil.create(getContext(), "test", 0, 1, BreakendDirection.Forward, B("TTT"), new byte[] { 1,2,3 }, B("ATTT"), new byte[] { 1,2,3,4}, 1, 1, 0);
		assertArrayEquals(new byte[] {1, 2, 3}, dba.getBreakpointQuality());
	}
	@Test
	public void should_set_assembly_evidence() {
		VariantContextDirectedBreakpoint dba = ReadEvidenceAssemblerUtil.create(getContext(), "test", 0, 1, BreakendDirection.Forward, B("TTT"), new byte[] { 1,2,3 }, B("ATTT"), new byte[] { 1,2,3,4}, 1, 2, 3);
		assertEquals(1, dba.getBreakendSummary().evidence.get(VcfAttributes.ASSEMBLY_READS));
		assertEquals(2, dba.getBreakendSummary().evidence.get(VcfAttributes.ASSEMBLY_BASES));
		assertEquals(4, dba.getBreakendSummary().evidence.get(VcfAttributes.ASSEMBLY_LENGTH));		
	}
	@Test
	public void should_set_assembly_attributes() {
		VariantContextDirectedBreakpoint dba = ReadEvidenceAssemblerUtil.create(getContext(), "test", 0, 1, BreakendDirection.Forward, B("TTT"), new byte[] { 1,2,3 }, B("ATTT"), new byte[] { 1,2,3,4}, 1, 2, 3);
		assertEquals("ATTT", dba.getAssemblyConsensus());
		assertEquals(3, dba.getAssemblyQuality(), 0);
		assertEquals("test", dba.getAssemblerProgram());
	}
}
