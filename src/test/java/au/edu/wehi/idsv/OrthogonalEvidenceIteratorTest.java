package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.Lists;

import au.edu.wehi.idsv.vcf.VcfFilter;


public class OrthogonalEvidenceIteratorTest extends TestHelper {
	@Test
	public void should_not_remove_unrelated_evidence() {
		DirectedEvidence e0 = SCE(FWD, Read(0, 1, "1M12S"));
		DirectedEvidence e1 = SCE(FWD, Read(0, 10, "1M10S"));
		DirectedEvidence e2 = SCE(FWD, Read(0, 10, "1M11S"));
		DirectedEvidence e3 = SCE(FWD, Read(0, 10, "1M12S"));
		DirectedEvidence e4 = SCE(FWD, Read(0, 100, "1M12S"));
		AssemblyEvidence a1 = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, Lists.transform(Lists.newArrayList(e1, e2, e3), EID), 0, 10, 1, B("TT"), B("TT"));
		List<DirectedEvidence> result = Lists.newArrayList(new OrthogonalEvidenceIterator(getContext().getLinear(), Lists.newArrayList(e0, e1, e2, e3, a1, e4).iterator(), 2, false));
		assertTrue(result.contains(e0));
		assertTrue(result.contains(e4));
		assertTrue(result.contains(a1));
	}
	@Test
	public void should_hydrate_all_assemblies() {
		DirectedEvidence e0 = SCE(FWD, Read(0, 1, "1M12S"));
		DirectedEvidence e1 = SCE(FWD, Read(0, 10, "1M10S"));
		DirectedEvidence e2 = SCE(FWD, Read(0, 10, "1M11S"));
		DirectedEvidence e3 = SCE(FWD, Read(0, 10, "1M12S"));
		DirectedEvidence e4 = SCE(FWD, Read(0, 100, "1M12S"));
		AssemblyEvidence a1 = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, Lists.transform(Lists.newArrayList(e2, e3), EID), 0, 10, 1, B("TT"), B("TT"));
		AssemblyEvidence a2 = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, Lists.transform(Lists.newArrayList(e1, e2), EID), 0, 11, 1, B("TT"), B("TT"));
		List<DirectedEvidence> result = Lists.newArrayList(new OrthogonalEvidenceIterator(getContext().getLinear(), Lists.newArrayList(e0, e1, e2, e3, a1, a2, e4).iterator(), 1000, true));
		assertEquals(2, result.size());
		assertEquals(2, ((SAMRecordAssemblyEvidence)result.get(0)).getEvidence().size());
		assertEquals(2, ((SAMRecordAssemblyEvidence)result.get(1)).getEvidence().size());
		assertTrue(((SAMRecordAssemblyEvidence)result.get(0)).getEvidence().contains(e3));
		assertTrue(((SAMRecordAssemblyEvidence)result.get(0)).getEvidence().contains(e2));
		assertTrue(((SAMRecordAssemblyEvidence)result.get(0)).getEvidence().contains(e3));
		assertTrue(((SAMRecordAssemblyEvidence)result.get(1)).getEvidence().contains(e1));
		assertTrue(((SAMRecordAssemblyEvidence)result.get(1)).getEvidence().contains(e2));
	}
	@Test
	public void assembliesOnly_should_not_output_nonassembly_evidence() {
		DirectedEvidence e0 = SCE(FWD, Read(0, 1, "1M12S"));
		DirectedEvidence e1 = SCE(FWD, Read(0, 10, "1M10S"));
		DirectedEvidence e2 = SCE(FWD, Read(0, 10, "1M11S"));
		DirectedEvidence e3 = SCE(FWD, Read(0, 10, "1M12S"));
		DirectedEvidence e4 = SCE(FWD, Read(0, 100, "1M12S"));
		AssemblyEvidence a1 = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, Lists.transform(Lists.newArrayList(e1, e2, e3), EID), 0, 10, 1, B("TT"), B("TT"));
		List<DirectedEvidence> result = Lists.newArrayList(new OrthogonalEvidenceIterator(getContext().getLinear(), Lists.newArrayList(e0, e1, e2, e3, a1, e4).iterator(), 2, true));
		assertEquals(1, result.size());
	}
	@Test
	public void should_remove_component_evidence() {
		DirectedEvidence e0 = SCE(FWD, Read(0, 1, "1M12S"));
		DirectedEvidence e1 = SCE(FWD, Read(0, 10, "1M10S"));
		DirectedEvidence e2 = SCE(FWD, Read(0, 10, "1M11S"));
		DirectedEvidence e3 = SCE(FWD, Read(0, 10, "1M12S"));
		DirectedEvidence e4 = SCE(FWD, Read(0, 100, "1M12S"));
		AssemblyEvidence a1 = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, Lists.transform(Lists.newArrayList(e1, e2, e3), EID), 0, 10, 1, B("TT"), B("TT"));
		List<DirectedEvidence> result = Lists.newArrayList(new OrthogonalEvidenceIterator(getContext().getLinear(), Lists.newArrayList(e0, e1, e2, e3, a1, e4).iterator(), 2, false));
		assertFalse(result.contains(e1));
		assertFalse(result.contains(e2));
		assertFalse(result.contains(e3));
	}
	@Test
	public void should_add_component_evidence_to_all_assemblies() {
		DirectedEvidence e0 = SCE(FWD, Read(0, 1, "1M12S"));
		DirectedEvidence e1 = SCE(FWD, Read(0, 10, "1M10S"));
		DirectedEvidence e2 = SCE(FWD, Read(0, 10, "1M11S"));
		DirectedEvidence e3 = SCE(FWD, Read(0, 10, "1M12S"));
		DirectedEvidence e4 = SCE(FWD, Read(0, 100, "1M12S"));
		SAMRecordAssemblyEvidence a1 = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, Lists.transform(Lists.newArrayList(e1, e2, e3), EID), 0, 10, 1, B("TT"), B("TT"));
		a1 = AssemblyFactory.hydrate(AES(), a1.getBackingRecord());
		List<DirectedEvidence> result = Lists.newArrayList(new OrthogonalEvidenceIterator(getContext().getLinear(), Lists.newArrayList(e0, e1, e2, e3, a1, e4).iterator(), 2, false));
		assertFalse(result.contains(e1));
		assertFalse(result.contains(e2));
		assertFalse(result.contains(e3));
		assertEquals(3, a1.getEvidence().size());
	}
	@Test
	public void should_not_remove_component_evidence_of_filtered_assembly() {
		DirectedEvidence e0 = SCE(FWD, Read(0, 1, "1M12S"));
		DirectedEvidence e1 = SCE(FWD, Read(0, 10, "1M10S"));
		DirectedEvidence e2 = SCE(FWD, Read(0, 10, "1M11S"));
		DirectedEvidence e3 = SCE(FWD, Read(0, 10, "1M12S"));
		DirectedEvidence e4 = SCE(FWD, Read(0, 100, "1M12S"));
		AssemblyEvidence a1 = AssemblyFactory.createAnchoredBreakend(getContext(), AES(), FWD, Lists.transform(Lists.newArrayList(e1, e2, e3), EID), 0, 10, 1, B("TT"), B("TT"));
		a1.filterAssembly(VcfFilter.SMALL_INDEL);
		List<DirectedEvidence> result = Lists.newArrayList(new OrthogonalEvidenceIterator(getContext().getLinear(), Lists.newArrayList(e0, e1, e2, e3, a1, e4).iterator(), 2, false));
		assertTrue(result.contains(e1));
		assertTrue(result.contains(e2));
		assertTrue(result.contains(e3));
	}
	@Test
	public void should_output_sorted() {
		ProcessingContext pc = getContext();
		SAMEvidenceSource ses = SES();
		AssemblyEvidenceSource aes = AES();
		List<DirectedEvidence> list = new ArrayList<DirectedEvidence>();
		for (int i = 1; i <= 1000; i++) {
			list.add(SoftClipEvidence.create(ses, FWD, Read(0, i, "1M10S"), null));
			list.add(AssemblyFactory.createAnchoredBreakend(pc, aes, FWD, null, 0, i, 1, B("TT"), B("TT")));
		}
		List<DirectedEvidence> result = Lists.newArrayList(new OrthogonalEvidenceIterator(pc.getLinear(), list.iterator(), 10, false));
		for (int i = 1; i < result.size(); i++) {
			assertTrue(DirectedEvidence.ByStartEnd.compare(result.get(i-1), result.get(i)) <= 0);
		}
	}
}
