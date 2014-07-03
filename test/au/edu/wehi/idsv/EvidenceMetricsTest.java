package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.fail;

import org.junit.Test;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.EvidenceMetrics;
import au.edu.wehi.idsv.NonReferenceReadPair;
import au.edu.wehi.idsv.ReadEvidenceAssemblerUtil;
import au.edu.wehi.idsv.SoftClipEvidence;
import au.edu.wehi.idsv.vcf.VcfAttributes;


public class EvidenceMetricsTest extends TestHelper {
	//@Test
	public void should_take_max_of_max_evidence() {
		fail();
	}
	@Test
	public void should_not_have_evidence_attribute_that_is_never_set() {
		EvidenceMetrics m = new EvidenceMetrics();
		// DP
		m.add(new NonReferenceReadPair(DP(0, 1, "100M", true, 1, 1, "100M", true)[0], DP(0, 1, "100M", true, 1, 1, "100M", true)[1], 300).getBreakendSummary().evidence);
		// OEA
		m.add(new NonReferenceReadPair(OEA(0, 1, "100M", true)[0], OEA(0, 1, "100M", true)[1], 300).getBreakendSummary().evidence);
		// SC
		m.add(new SoftClipEvidence(new SoftClipEvidence(getContext(), BreakendDirection.Backward, Read(0, 1, "6S4M")), withMapq(1,Read(0, 1, "1S3M2S"))[0]).getBreakendSummary().evidence);
		// Assembly
		m.add(ReadEvidenceAssemblerUtil.breakendBuilder(getContext(), "test", 0, 10, BreakendDirection.Backward, B("AT"), B("ATA"), 1, 1, 1).make().getBreakendSummary().evidence);
		for (VcfAttributes a : VcfAttributes.evidenceValues()) {
			assertFalse(a.name(), m.get(a) == 0);
		}
	}
	@Test
	public void add_should_add_evidence() {
		EvidenceMetrics e = new EvidenceMetrics();
		e.set(VcfAttributes.ASSEMBLY_LENGTH, 5);
		EvidenceMetrics a = new EvidenceMetrics();
		a.set(VcfAttributes.SOFT_CLIP_READ_COUNT, 2);
		a.set(VcfAttributes.ASSEMBLY_BASES, 1);
		e.add(a);
		e.add(a);
		assertEquals(4, e.get(VcfAttributes.SOFT_CLIP_READ_COUNT));
		assertEquals(2, e.get(VcfAttributes.ASSEMBLY_BASES));
		assertEquals(5, e.get(VcfAttributes.ASSEMBLY_LENGTH));
	}
	@Test
	public void remove_should_subtract_evidence() {
		EvidenceMetrics e = new EvidenceMetrics();
		e.set(VcfAttributes.SOFT_CLIP_READ_COUNT, 5);
		e.set(VcfAttributes.ASSEMBLY_BASES, 7);
		e.set(VcfAttributes.ASSEMBLY_LENGTH, 9);
		EvidenceMetrics a = new EvidenceMetrics();
		a.set(VcfAttributes.SOFT_CLIP_READ_COUNT, 2);
		a.set(VcfAttributes.ASSEMBLY_BASES, 1);
		e.remove(a);
		e.remove(a);
		assertEquals(1, e.get(VcfAttributes.SOFT_CLIP_READ_COUNT));
		assertEquals(5, e.get(VcfAttributes.ASSEMBLY_BASES));
		assertEquals(9, e.get(VcfAttributes.ASSEMBLY_LENGTH));
	}
	@Test
	public void constructor_should_deep_copy() {
		EvidenceMetrics e = new EvidenceMetrics();
		e.set(VcfAttributes.SOFT_CLIP_READ_COUNT, 5);
		EvidenceMetrics a = new EvidenceMetrics(e);
		e.set(VcfAttributes.SOFT_CLIP_READ_COUNT, 2);
		a.set(VcfAttributes.SOFT_CLIP_READ_COUNT, 3);
		
		assertEquals(2, e.get(VcfAttributes.SOFT_CLIP_READ_COUNT));
		assertEquals(3, a.get(VcfAttributes.SOFT_CLIP_READ_COUNT));
	}
	@Test
	public void clone_should_deep_copy() {
		EvidenceMetrics e = new EvidenceMetrics();
		e.set(VcfAttributes.SOFT_CLIP_READ_COUNT, 5);
		EvidenceMetrics a = e.clone();
		e.set(VcfAttributes.SOFT_CLIP_READ_COUNT, 2);
		a.set(VcfAttributes.SOFT_CLIP_READ_COUNT, 3);
		
		assertEquals(2, e.get(VcfAttributes.SOFT_CLIP_READ_COUNT));
		assertEquals(3, a.get(VcfAttributes.SOFT_CLIP_READ_COUNT));
	}
}
