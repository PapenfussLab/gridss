package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import com.google.common.collect.ImmutableList;

import au.edu.wehi.idsv.sam.SamTags;
import htsjdk.samtools.SAMRecord;

public class AssemblyAttributesTest extends TestHelper {
	@Test
	public void setEvidenceIDs_should_set_same_value_regardless_of_collection_ordering() {
		DirectedEvidence e1 = SCE(FWD, withMapq(10, Read(0, 1, "1M2S")));
		DirectedEvidence e2 = SCE(FWD, withMapq(20, Read(0, 1, "1M1S")));
		SAMRecord ass1 = AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(0, FWD, 1, 1, 2), ImmutableList.of(e1, e2), B("GTAC"), new byte[] {1,2,3,4}, new int[] {0, 0});
		SAMRecord ass2 = AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(0, FWD, 1, 1, 2), ImmutableList.of(e2, e1), B("GTAC"), new byte[] {1,2,3,4}, new int[] {0, 0});
		
		assertEquals(ass1.getAttribute(SamTags.EVIDENCEID), ass2.getAttribute(SamTags.EVIDENCEID));
	}
	@Test
	public void isPartOfAssembly_should_use_evidenceId() {
		DirectedEvidence e1 = SCE(FWD, withMapq(10, Read(0, 1, "1M2S")));
		DirectedEvidence e2 = SCE(FWD, withMapq(20, Read(0, 2, "1M1S")));
		DirectedEvidence e3 = SCE(FWD, withMapq(20, Read(0, 3, "1M1S")));
		SAMRecord ass1 = AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), new SequentialIdGenerator("asm"), new BreakendSummary(0, FWD, 1, 1, 2), ImmutableList.of(e1, e2), B("GTAC"), new byte[] {1,2,3,4}, new int[] {0, 0});
		
		AssemblyAttributes attr = new AssemblyAttributes(ass1);
		assertTrue(attr.isPartOfAssembly(e1));
		assertTrue(attr.isPartOfAssembly(e2));
		assertFalse(attr.isPartOfAssembly(e3));
	}
}
