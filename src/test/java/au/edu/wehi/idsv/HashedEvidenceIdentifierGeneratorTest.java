package au.edu.wehi.idsv;

import java.util.regex.Pattern;

import org.junit.Assert;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class HashedEvidenceIdentifierGeneratorTest extends EvidenceIdentifierGeneratorTest {
	@Override
	public EvidenceIdentifierGenerator getGeneratorToTest() {
		return new HashedEvidenceIdentifierGenerator();
	}
	@Test
	public void should_default_to_32_byte_hash() {
		SAMRecord r = Read(0, 1, "10M10S");
		SoftClipEvidence sc = SCE(FWD, r);
		HashedEvidenceIdentifierGenerator gen = new HashedEvidenceIdentifierGenerator();
		Assert.assertEquals(32, gen.getEvidenceID(sc).length());
	}
	@Test
	public void hash_should_use_Base64_url_encoding() {
		HashedEvidenceIdentifierGenerator gen = new HashedEvidenceIdentifierGenerator();
		for (int i = 0; i < 1024; i++) {
			SAMRecord r = Read(0, 1, "10M10S");
			r.setReadName(String.format("read%d", i));
			SoftClipEvidence sc = SCE(FWD, r);
			String hash = gen.getEvidenceID(sc);
			Assert.assertTrue(Pattern.matches("^[0-9a-zA-Z_-]{32}$", hash));
		}
	}
}
