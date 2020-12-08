package au.edu.wehi.idsv;

import au.edu.wehi.idsv.sam.SamTags;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Range;
import htsjdk.samtools.SAMRecord;
import org.junit.Test;

import static org.junit.Assert.*;

public class AssemblyAttributesTest extends TestHelper {
	@Test
	public void setEvidenceIDs_should_set_same_value_regardless_of_collection_ordering() {
		DirectedEvidence e1 = SCE(FWD, withMapq(10, Read(0, 1, "1M2S")));
		DirectedEvidence e2 = SCE(FWD, withMapq(20, Read(0, 1, "1M1S")));
		SAMRecord ass1 = AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(0, FWD, 1, 1, 2), ImmutableList.of(e1, e2), null, B("GTAC"), new byte[] {1,2,3,4});
		SAMRecord ass2 = AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(0, FWD, 1, 1, 2), ImmutableList.of(e2, e1), null, B("GTAC"), new byte[] {1,2,3,4});
		
		assertEquals(ass1.getAttribute(SamTags.ASSEMBLY_EVIDENCE_EVIDENCEID ), ass2.getAttribute(SamTags.ASSEMBLY_EVIDENCE_EVIDENCEID));
	}
	@Test
	public void isPartOfAssembly_should_use_evidenceId() {
		DirectedEvidence e1 = SCE(FWD, withMapq(10, Read(0, 1, "1M2S")));
		DirectedEvidence e2 = SCE(FWD, withMapq(20, Read(0, 2, "1M1S")));
		DirectedEvidence e3 = SCE(FWD, withMapq(20, Read(0, 3, "1M1S")));
		SAMRecord ass1 = AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(0, FWD, 1, 1, 2), ImmutableList.of(e1, e2), null, B("GTAC"), new byte[] {1,2,3,4});
		
		AssemblyAttributes attr = new AssemblyAttributes(ass1);
		assertTrue(attr.isPartOfAssembly(e1));
		assertTrue(attr.isPartOfAssembly(e2));
		assertFalse(attr.isPartOfAssembly(e3));
	}

	@Test
	public void getEvidenceIDs_should_subset_on_assembly_offset() {
		DirectedEvidence e1 = SCE(FWD, withMapq(10, Read(0, 1, "1M2S")));
		DirectedEvidence e2 = SCE(FWD, withMapq(20, Read(0, 2, "1M1S")));
		DirectedEvidence e3 = SCE(FWD, withMapq(20, Read(0, 3, "1M1S")));
		SAMRecord ass1 = AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(0, FWD, 1, 1, 2), ImmutableList.of(e1, e2), null, B("GTAC"), new byte[] {1,2,3,4});

		AssemblyAttributes attr = new AssemblyAttributes(ass1);
		assertTrue(attr.isPartOfAssembly(e1));
		assertTrue(attr.isPartOfAssembly(e2));
		assertFalse(attr.isPartOfAssembly(e3));
	}
	@Test
	public void getMinQualPosition_should_return_min_position() {
		SAMRecord r = Read(0, 1, "100M");
		r.setAttribute(SamTags.IS_ASSEMBLY, 1);
		r.setAttribute(SamTags.ASSEMBLY_DIRECTION, "f");
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_TYPE, new byte[] { 0,0, 0});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_CATEGORY, new int[] { 0, 0, 0});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_OFFSET_START, new int[] { 1, 2, 3});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_OFFSET_END, new int[] { 10, 11, 5});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_QUAL, new float[] { 1, 2, 4});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_EVIDENCEID, "1 2 3");
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_FRAGMENTID, "1 2 3");
		AssemblyAttributes attr = new AssemblyAttributes(r);
		assertEquals(12, attr.getMinQualPosition(Range.closed(1, 20), null, null, null));
	}
	@Test
	public void getMaxQualPosition_should_return_min_position() {
		SAMRecord r = Read(0, 1, "100M");
		r.setAttribute(SamTags.IS_ASSEMBLY, 1);
		r.setAttribute(SamTags.ASSEMBLY_DIRECTION, "f");
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_TYPE, new byte[] { 0,0, 0});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_CATEGORY, new int[] { 0, 0, 0});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_OFFSET_START, new int[] { 1, 2, 3});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_OFFSET_END, new int[] { 10, 11, 5});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_QUAL, new float[] { 1, 2, 4});
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_EVIDENCEID, "1 2 3");
		r.setAttribute(SamTags.ASSEMBLY_EVIDENCE_FRAGMENTID, "1 2 3");
		AssemblyAttributes attr = new AssemblyAttributes(r);
		assertEquals(3, attr.getMaxQualPosition(Range.closed(0, 15), null, null, null));
	}
}
