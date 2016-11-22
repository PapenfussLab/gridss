package au.edu.wehi.idsv;

import static org.junit.Assert.*;

import java.util.List;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class SingleReadEvidenceTest extends TestHelper {
	@Test
	public void should_not_create_indel_for_XNX_placeholder() {
		for (SAMRecord r : new SAMRecord[] {
				Read(0, 1, "1X10S"),
				Read(0, 1, "2X10S"),
				Read(0, 1, "1X1N1X10S"),
				Read(0, 1, "1X100N1X10S"),
		}) {
			List<SingleReadEvidence> list = SingleReadEvidence.createEvidence(SES(), r);
			assertEquals(1, list.size());
			assertTrue(list.get(0) instanceof SoftClipEvidence);
		}
	}
}
