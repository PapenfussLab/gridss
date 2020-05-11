package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class StringEvidenceIdentifierGeneratorTest extends EvidenceIdentifierGeneratorTest {
	@Override
	public EvidenceIdentifierGenerator getGeneratorToTest() {
		return new StringEvidenceIdentifierGenerator();
	}
	@Test
	public void should_use_hash_separator() {
		SAMRecord r = withName("readname", Read(0, 1, "5M1D1M4S"))[0];
		SoftClipEvidence sce = SCE(FWD, ses, r);
		assertEquals("readname#0", gen.extractSegmentUniqueName(sce.getEvidenceID()));
		assertEquals("readname#0#polyA#1#+#5M1D1M4S", gen.extractAlignmentUniqueName(sce.getEvidenceID()));
	}
}
