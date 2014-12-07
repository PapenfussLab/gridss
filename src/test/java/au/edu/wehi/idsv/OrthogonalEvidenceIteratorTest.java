package au.edu.wehi.idsv;

import static org.junit.Assert.*;

import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.vcf.VcfFilter;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;


public class OrthogonalEvidenceIteratorTest extends TestHelper {
	@Test
	public void should_not_remove_unrelated_evidence() {
		DirectedEvidence e0 = SCE(FWD, Read(0, 1, "1M12S"));
		DirectedEvidence e1 = SCE(FWD, Read(0, 10, "1M10S"));
		DirectedEvidence e2 = SCE(FWD, Read(0, 10, "1M11S"));
		DirectedEvidence e3 = SCE(FWD, Read(0, 10, "1M12S"));
		DirectedEvidence e4 = SCE(FWD, Read(0, 100, "1M12S"));
		AssemblyEvidence a1 = AssemblyFactory.createAnchored(getContext(), AES(), FWD, Sets.newHashSet(e1, e2, e3), 0, 10, 1, B("TT"), B("TT"), 1, 1);
		List<DirectedEvidence> result = Lists.newArrayList(new OrthogonalEvidenceIterator(getContext().getLinear(), Lists.newArrayList(e0, e1, e2, e3, a1, e4).iterator(), 2));
		assertTrue(result.contains(e0));
		assertTrue(result.contains(e4));
		assertTrue(result.contains(a1));
	}
	@Test
	public void should_remove_component_evidence() {
		DirectedEvidence e0 = SCE(FWD, Read(0, 1, "1M12S"));
		DirectedEvidence e1 = SCE(FWD, Read(0, 10, "1M10S"));
		DirectedEvidence e2 = SCE(FWD, Read(0, 10, "1M11S"));
		DirectedEvidence e3 = SCE(FWD, Read(0, 10, "1M12S"));
		DirectedEvidence e4 = SCE(FWD, Read(0, 100, "1M12S"));
		AssemblyEvidence a1 = AssemblyFactory.createAnchored(getContext(), AES(), FWD, Sets.newHashSet(e1, e2, e3), 0, 10, 1, B("TT"), B("TT"), 1, 1);
		List<DirectedEvidence> result = Lists.newArrayList(new OrthogonalEvidenceIterator(getContext().getLinear(), Lists.newArrayList(e0, e1, e2, e3, a1, e4).iterator(), 2));
		assertFalse(result.contains(e1));
		assertFalse(result.contains(e2));
		assertFalse(result.contains(e3));
	}
	@Test
	public void should_not_remove_component_evidence_of_filtered_assembly() {
		DirectedEvidence e0 = SCE(FWD, Read(0, 1, "1M12S"));
		DirectedEvidence e1 = SCE(FWD, Read(0, 10, "1M10S"));
		DirectedEvidence e2 = SCE(FWD, Read(0, 10, "1M11S"));
		DirectedEvidence e3 = SCE(FWD, Read(0, 10, "1M12S"));
		DirectedEvidence e4 = SCE(FWD, Read(0, 100, "1M12S"));
		AssemblyEvidence a1 = AssemblyFactory.createAnchored(getContext(), AES(), FWD, Sets.newHashSet(e1, e2, e3), 0, 10, 1, B("TT"), B("TT"), 1, 1);
		a1.filterAssembly(VcfFilter.SMALL_INDEL);
		List<DirectedEvidence> result = Lists.newArrayList(new OrthogonalEvidenceIterator(getContext().getLinear(), Lists.newArrayList(e0, e1, e2, e3, a1, e4).iterator(), 2));
		assertTrue(result.contains(e1));
		assertTrue(result.contains(e2));
		assertTrue(result.contains(e3));
	}
}
