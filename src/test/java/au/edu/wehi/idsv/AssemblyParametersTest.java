package au.edu.wehi.idsv;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import au.edu.wehi.idsv.vcf.VcfFilter;

import com.google.common.collect.Sets;


public class AssemblyParametersTest extends TestHelper {
	@Test
	public void should_filter_fully_reference_assemblies() {
		AssemblyEvidence e = AssemblyFactory.createAnchoredBreakend(
				getContext(), AES(), BWD, Sets.<DirectedEvidence>newHashSet(),
				0, 1, 2, B("AA"), B("AA"), 2, 0);
		getContext().getAssemblyParameters().applyFilters(e);
		assertTrue(e.isAssemblyFiltered());
		assertTrue(e.getFilters().contains(VcfFilter.REFERENCE_ALLELE));
	}
	@Test
	public void should_filter_single_read_assemblies() {
		AssemblyEvidence e = AssemblyFactory.createAnchoredBreakend(
				getContext(), AES(), FWD, Sets.<DirectedEvidence>newHashSet(
					SCE(FWD, Read(0, 1, "1M2S"))),
				0, 1, 1, B("AA"), B("AA"), 2, 0);
		getContext().getAssemblyParameters().applyFilters(e);
		assertTrue(e.isAssemblyFiltered());
		assertTrue(e.getFilters().contains(VcfFilter.ASSEMBLY_TOO_FEW_READ));
	}
	@Test
	public void should_filter_mate_anchored_assembly_shorter_than_read_length() {
		AssemblyEvidence e = AssemblyFactory.createUnanchoredBreakend(
				getContext(), AES(), Sets.<DirectedEvidence>newHashSet(NRRP(OEA(0, 1, "4M", true))),
				B("AAA"), B("AAA"), 2, 0);
		getContext().getAssemblyParameters().applyFilters(e);
		assertTrue(e.isAssemblyFiltered());
		assertTrue(e.getFilters().contains(VcfFilter.ASSEMBLY_TOO_SHORT));
	}
	@Test
	public void read_length_filter_should_not_apply_to_anchored_breakend_assembly() {
		AssemblyEvidence e = AssemblyFactory.createAnchoredBreakend(
				getContext(), AES(), FWD, Sets.<DirectedEvidence>newHashSet(
						SCE(FWD, Read(0, 1, "1M2S")),
						NRRP(OEA(0, 1, "3M", true)),
						NRRP(OEA(0, 1, "4M", true))),
				0, 1, 1, B("AAA"), B("AAA"), 2, 0);
		getContext().getAssemblyParameters().applyFilters(e);
		assertTrue(!e.getFilters().contains(VcfFilter.ASSEMBLY_TOO_SHORT));
	}
	@Test
	public void soft_clip_size_filter_should_not_apply_to_unanchored_assembly() {
		AssemblyEvidence e = AssemblyFactory.createUnanchoredBreakend(
				getContext(), AES(), Sets.<DirectedEvidence>newHashSet(
						NRRP(OEA(0, 1, "3M", false)),
						NRRP(OEA(0, 1, "4M", false))),
				B("AAAAAA"), B("AAAAAA"), 2, 0);
		getContext().getAssemblyParameters().applyFilters(e);
		assertFalse(e.isAssemblyFiltered());
	}
	@Test
	public void should_filter_if_no_breakpoint_assembly() {
		AssemblyEvidence e = AssemblyFactory.createAnchoredBreakend(
				getContext(), AES(), BWD, Sets.<DirectedEvidence>newHashSet(
						SCE(BWD, Read(0, 1, "1S1M")),
						NRRP(OEA(0, 1, "3M", false)),
						NRRP(OEA(0, 1, "4M", false))),
				0, 1, 2, B("AA"), B("AA"), 2, 0);
		assertTrue(getContext().getAssemblyParameters().applyFilters(e));
		assertTrue(e.getFilters().contains(VcfFilter.REFERENCE_ALLELE));
	}
	@Test
	public void should_not_apply_breakend_filter_to_unanchored_assembly() {
		AssemblyEvidence e = AssemblyFactory.createUnanchoredBreakend(
				getContext(), AES(), Sets.<DirectedEvidence>newHashSet(
						NRRP(SES(100, 100), DP(0, 1, "1M", true, 0, 5, "1M", false))),
				B("AA"), B("AA"), 2, 0);
		getContext().getAssemblyParameters().applyFilters(e);
		assertFalse(e.getFilters().contains(VcfFilter.REFERENCE_ALLELE));
	}
	@Test
	public void should_filter_too_few_reads() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().minReads = 3;
		AssemblyEvidence e = AssemblyFactory.createAnchoredBreakend(pc, AES(), BWD, Sets.<DirectedEvidence>newHashSet(), 0, 1, 5, B("AACGTG"), B("AACGTG"), 2, 0);
		assertTrue(pc.getAssemblyParameters().applyFilters(e));
		assertTrue(e.getFilters().contains(VcfFilter.ASSEMBLY_TOO_FEW_READ));
		
		e = AssemblyFactory.createAnchoredBreakend(
				getContext(), AES(), BWD, Sets.<DirectedEvidence>newHashSet(
						SCE(BreakendDirection.Forward, withQual(new byte[] { 5,5,5,5,5,5 }, withSequence("AACGTG", Read(0, 1, "1M5S"))))),
				0, 1, 5, B("AACGTG"), B("AACGTG"), 2, 0);
		assertTrue(pc.getAssemblyParameters().applyFilters(e));
		assertTrue(e.getFilters().contains(VcfFilter.ASSEMBLY_TOO_FEW_READ));
		
		e = AssemblyFactory.createAnchoredBreakend(
				pc, AES(), BWD, Sets.<DirectedEvidence>newHashSet(
						SCE(BreakendDirection.Forward, withQual(new byte[] { 5,5,5,5,5,5 }, withSequence("AACGTG", Read(0, 1, "1M5S")))),
						NRRP(OEA(0, 1, "4M", false))),
				0, 1, 5, B("AACGTG"), B("AACGTG"), 2, 0);
		assertTrue(pc.getAssemblyParameters().applyFilters(e));
		assertTrue(e.getFilters().contains(VcfFilter.ASSEMBLY_TOO_FEW_READ));
		
		e = AssemblyFactory.createAnchoredBreakend(
				pc, AES(), BWD, Sets.<DirectedEvidence>newHashSet(
						SCE(BreakendDirection.Forward, withQual(new byte[] { 5,5,5,5,5,5 }, withSequence("AACGTG", Read(0, 1, "1M5S")))),
						NRRP(OEA(0, 1, "3M", false)),
						NRRP(OEA(0, 1, "5M", false))),
				0, 1, 5, B("AACGTG"), B("AACGTG"), 2, 0);
		assertFalse(e.getFilters().contains(VcfFilter.ASSEMBLY_TOO_FEW_READ));
	}
}
