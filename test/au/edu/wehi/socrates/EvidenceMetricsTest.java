package au.edu.wehi.socrates;

import static org.junit.Assert.*;
import htsjdk.samtools.SAMRecord;

import org.junit.Test;

import au.edu.wehi.socrates.vcf.EvidenceAttributes;


public class EvidenceMetricsTest extends TestHelper {
	//@Test
	public void should_take_max_of_max_evidence() {
		fail();
	}
	@Test
	public void should_not_have_evidence_attribute_that_is_never_set() {
		EvidenceMetrics m = new EvidenceMetrics();
		// DP
		m.add(new NonReferenceReadPair(DP(0, 1, "100M", true, 0, 1, "100M", true)[0], DP(0, 1, "100M", true, 0, 1, "100M", true)[1], 100).getBreakendSummary().evidence);
		// OEA
		m.add(new NonReferenceReadPair(OEA(0, 1, "100M", true)[0], OEA(0, 1, "100M", true)[1], 100).getBreakendSummary().evidence);
		// SC
		m.add(new SoftClipEvidence(new SoftClipEvidence(getContext(), BreakendDirection.Backward, Read(0, 1, "6S4M")), withMapq(1,Read(0, 1, "1S3M2S"))[0]).getBreakendSummary().evidence);
		// Assembly
		m.add(ReadEvidenceAssemblerUtil.create(getContext(), "test", 0, 10, BreakendDirection.Backward, B("AT"), B("ATA"), 1, 1, 1).getBreakendSummary().evidence);
		for (EvidenceAttributes a : EvidenceAttributes.values()) {
			assertNotEquals(a.name(), 0, m.get(a));
		}
	}
	@Test
	public void add_should_add_evidence() {
		EvidenceMetrics e = new EvidenceMetrics();
		e.set(EvidenceAttributes.ASSEMBLY_LENGTH, 5);
		EvidenceMetrics a = new EvidenceMetrics();
		a.set(EvidenceAttributes.SOFT_CLIP_READ_COUNT, 2);
		a.set(EvidenceAttributes.ASSEMBLY_BASES, 1);
		e.add(a);
		e.add(a);
		assertEquals(4, e.get(EvidenceAttributes.SOFT_CLIP_READ_COUNT));
		assertEquals(2, e.get(EvidenceAttributes.ASSEMBLY_BASES));
		assertEquals(5, e.get(EvidenceAttributes.ASSEMBLY_LENGTH));
	}
	@Test
	public void remove_should_subtract_evidence() {
		EvidenceMetrics e = new EvidenceMetrics();
		e.set(EvidenceAttributes.SOFT_CLIP_READ_COUNT, 5);
		e.set(EvidenceAttributes.ASSEMBLY_BASES, 7);
		e.set(EvidenceAttributes.ASSEMBLY_LENGTH, 9);
		EvidenceMetrics a = new EvidenceMetrics();
		a.set(EvidenceAttributes.SOFT_CLIP_READ_COUNT, 2);
		a.set(EvidenceAttributes.ASSEMBLY_BASES, 1);
		e.remove(a);
		e.remove(a);
		assertEquals(1, e.get(EvidenceAttributes.SOFT_CLIP_READ_COUNT));
		assertEquals(5, e.get(EvidenceAttributes.ASSEMBLY_BASES));
		assertEquals(9, e.get(EvidenceAttributes.ASSEMBLY_LENGTH));
	}
	@Test
	public void constructor_should_deep_copy() {
		EvidenceMetrics e = new EvidenceMetrics();
		e.set(EvidenceAttributes.SOFT_CLIP_READ_COUNT, 5);
		EvidenceMetrics a = new EvidenceMetrics(e);
		e.set(EvidenceAttributes.SOFT_CLIP_READ_COUNT, 2);
		a.set(EvidenceAttributes.SOFT_CLIP_READ_COUNT, 3);
		
		assertEquals(2, e.get(EvidenceAttributes.SOFT_CLIP_READ_COUNT));
		assertEquals(3, a.get(EvidenceAttributes.SOFT_CLIP_READ_COUNT));
	}
	@Test
	public void clone_should_deep_copy() {
		EvidenceMetrics e = new EvidenceMetrics();
		e.set(EvidenceAttributes.SOFT_CLIP_READ_COUNT, 5);
		EvidenceMetrics a = e.clone();
		e.set(EvidenceAttributes.SOFT_CLIP_READ_COUNT, 2);
		a.set(EvidenceAttributes.SOFT_CLIP_READ_COUNT, 3);
		
		assertEquals(2, e.get(EvidenceAttributes.SOFT_CLIP_READ_COUNT));
		assertEquals(3, a.get(EvidenceAttributes.SOFT_CLIP_READ_COUNT));
	}
}
