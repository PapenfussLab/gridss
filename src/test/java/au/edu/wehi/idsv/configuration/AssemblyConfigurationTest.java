package au.edu.wehi.idsv.configuration;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;

import org.junit.Test;

import au.edu.wehi.idsv.AssemblyFactory;
import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakendSummary;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.MockDirectedEvidence;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.vcf.VcfFilter;

import com.google.common.collect.Lists;


public class AssemblyConfigurationTest extends TestHelper {
	@Test
	public void should_filter_fully_reference_assemblies() {
		SAMRecordAssemblyEvidence e = AssemblyFactory.createAnchoredBreakend(
				getContext(), AES(), BWD, null,
				0, 1, 2, B("AA"), B("AA"))
				.annotateAssembly();
		getContext().getAssemblyParameters().applyAnnotationFilters(e);
		assertTrue(e.isAssemblyFiltered());
		assertTrue(e.getFilters().contains(VcfFilter.REFERENCE_ALLELE));
	}
	@Test
	public void should_filter_single_read_assemblies() {
		MockDirectedEvidence ev = new MockDirectedEvidence(new BreakendSummary(0, FWD, 1, 2));
		SAMRecordAssemblyEvidence e = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, Lists.newArrayList(ev.getEvidenceID()), 0, 1, 1, B("AA"), B("AA"))
				.hydrateEvidenceSet(ev)
				.annotateAssembly();
		getContext().getAssemblyParameters().applyAnnotationFilters(e);
		assertTrue(e.isAssemblyFiltered());
		assertTrue(e.getFilters().contains(VcfFilter.ASSEMBLY_TOO_FEW_READ));
	}
	@Test
	public void should_filter_mate_anchored_assembly_shorter_than_read_length() {
		ArrayList<DirectedEvidence> support = Lists.<DirectedEvidence>newArrayList(NRRP(OEA(0, 1, "4M", true)));
		SAMRecordAssemblyEvidence e = AssemblyFactory.createUnanchoredBreakend(
				getContext(), AES(), new BreakendSummary(0, FWD, 5, 300), Lists.transform(support, EID),
				B("AAA"), B("AAA"), new int[] { 2, 0 });
		e.hydrateEvidenceSet(support);
		e.annotateAssembly();
		getContext().getAssemblyParameters().applyAnnotationFilters(e);
		assertTrue(e.isAssemblyFiltered());
		assertTrue(e.getFilters().contains(VcfFilter.ASSEMBLY_TOO_SHORT));
	}
	@Test
	public void read_length_filter_should_not_apply_to_anchored_breakend_assembly() {
		ArrayList<DirectedEvidence> support = Lists.<DirectedEvidence>newArrayList(SCE(FWD, Read(0, 1, "1M2S")),
				NRRP(OEA(0, 1, "3M", true)),
				NRRP(OEA(0, 1, "4M", true)));
		SAMRecordAssemblyEvidence e = AssemblyFactory.createAnchoredBreakend(
				getContext(), AES(), FWD, Lists.transform(support, EID),
				0, 1, 1, B("AAA"), B("AAA"));
		e.hydrateEvidenceSet(support);
		e.annotateAssembly();
		getContext().getAssemblyParameters().applyAnnotationFilters(e);
		assertTrue(!e.getFilters().contains(VcfFilter.ASSEMBLY_TOO_SHORT));
	}
	@Test
	public void soft_clip_size_filter_should_not_apply_to_unanchored_assembly() {
		ArrayList<DirectedEvidence> support = Lists.<DirectedEvidence>newArrayList(
				NRRP(OEA(0, 1, "3M", false)),
				NRRP(OEA(0, 1, "4M", false)));
		SAMRecordAssemblyEvidence e = AssemblyFactory.createUnanchoredBreakend(
				getContext(), AES(), new BreakendSummary(0, FWD, 5, 300), Lists.transform(support, EID),
				B("AAAAAA"), B("AAAAAA"), new int[] { 2, 0 });
		e.hydrateEvidenceSet(support);
		e.annotateAssembly();
		getContext().getAssemblyParameters().applyAnnotationFilters(e);
		assertFalse(e.isAssemblyFiltered());
	}
	@Test
	public void should_filter_if_no_breakpoint_assembly() {
		ArrayList<DirectedEvidence> support = Lists.<DirectedEvidence>newArrayList(
				SCE(BWD, Read(0, 1, "1S1M")),
				NRRP(OEA(0, 1, "3M", false)),
				NRRP(OEA(0, 1, "4M", false)));
		SAMRecordAssemblyEvidence e = AssemblyFactory.createAnchoredBreakend(
				getContext(), AES(), BWD, Lists.transform(support, EID),
				0, 1, 2, B("AA"), B("AA"));
		e.hydrateEvidenceSet(support);
		e.annotateAssembly();
		assertTrue(getContext().getAssemblyParameters().applyAnnotationFilters(e));
		assertTrue(e.getFilters().contains(VcfFilter.REFERENCE_ALLELE));
	}
	@Test
	public void should_not_apply_breakend_filter_to_unanchored_assembly() {
		ArrayList<DirectedEvidence> support = Lists.<DirectedEvidence>newArrayList(NRRP(SES(100, 100), DP(0, 1, "1M", true, 0, 5, "1M", false)));
		SAMRecordAssemblyEvidence e = AssemblyFactory.createUnanchoredBreakend(
				getContext(), AES(), new BreakendSummary(0, FWD, 1, 300), Lists.transform(support, EID),
				B("AA"), B("AA"), new int[] { 2, 0});
		e.hydrateEvidenceSet(support);
		e.annotateAssembly();
		getContext().getAssemblyParameters().applyAnnotationFilters(e);
		assertFalse(e.getFilters().contains(VcfFilter.REFERENCE_ALLELE));
	}
	@Test
	public void should_filter_too_few_reads() {
		ProcessingContext pc = getContext();
		pc.getAssemblyParameters().minReads = 3;
		SAMRecordAssemblyEvidence e = AssemblyFactory.createAnchoredBreakend(pc, AES(), BWD, null, 0, 1, 5, B("AACGTG"), B("AACGTG"));
		e.annotateAssembly();
		assertTrue(pc.getAssemblyParameters().applyAnnotationFilters(e));
		assertTrue(e.getFilters().contains(VcfFilter.ASSEMBLY_TOO_FEW_READ));
		
		ArrayList<DirectedEvidence> support = Lists.<DirectedEvidence>newArrayList(
				SCE(BreakendDirection.Forward, withQual(new byte[] { 5,5,5,5,5,5 }, withSequence("AACGTG", Read(0, 1, "1M5S")))));
		e = AssemblyFactory.createAnchoredBreakend(
				getContext(), AES(), BWD, Lists.transform(support, EID),
				0, 1, 5, B("AACGTG"), B("AACGTG"));
		e.hydrateEvidenceSet(support);
		e.annotateAssembly();
		assertTrue(pc.getAssemblyParameters().applyAnnotationFilters(e));
		assertTrue(e.getFilters().contains(VcfFilter.ASSEMBLY_TOO_FEW_READ));
		
		support = Lists.<DirectedEvidence>newArrayList(
				SCE(BreakendDirection.Forward, withQual(new byte[] { 5,5,5,5,5,5 }, withSequence("AACGTG", Read(0, 1, "1M5S")))),
				NRRP(OEA(0, 1, "4M", false)));
		e = AssemblyFactory.createAnchoredBreakend(
				pc, AES(), BWD, Lists.transform(support, EID),
				0, 1, 5, B("AACGTG"), B("AACGTG"));
		e.hydrateEvidenceSet(support);
		e.annotateAssembly();
		assertTrue(pc.getAssemblyParameters().applyAnnotationFilters(e));
		assertTrue(e.getFilters().contains(VcfFilter.ASSEMBLY_TOO_FEW_READ));
		
		support = Lists.<DirectedEvidence>newArrayList(
				SCE(BreakendDirection.Forward, withQual(new byte[] { 5,5,5,5,5,5 }, withSequence("AACGTG", Read(0, 1, "1M5S")))),
				NRRP(OEA(0, 1, "3M", false)),
				NRRP(OEA(0, 1, "5M", false)));
		e = AssemblyFactory.createAnchoredBreakend(
				pc, AES(), BWD, Lists.transform(support, EID),
				0, 1, 5, B("AACGTG"), B("AACGTG"));
		e.hydrateEvidenceSet(support);
		e.annotateAssembly();
		assertFalse(e.getFilters().contains(VcfFilter.ASSEMBLY_TOO_FEW_READ));
	}
}
