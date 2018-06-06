package au.edu.wehi.idsv;

import java.util.List;

import org.junit.Assert;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;

public class SplitReadFastqExtractorTest extends TestHelper {
	@Test
	public void should_extract_full_sequence() {
		SplitReadFastqExtractor srfe = new SplitReadFastqExtractor(false, 1, 0, false, false, true, new StringEvidenceIdentifierGenerator());
		List<FastqRecord> result = srfe.extract(Read(0, 1, "25M75S"));
		Assert.assertEquals(1, result.size());
		Assert.assertEquals(100, result.get(0).getReadLength());
	}
	@Test
	public void should_extract_remaining_sequence() {
		SplitReadFastqExtractor srfe = new SplitReadFastqExtractor(false, 1, 0, false, false, true, new StringEvidenceIdentifierGenerator());
		FastqRecord fqr = srfe.extract(Read(0, 1, "25M75S")).get(0);
		SAMRecord aligned = new SAMRecord(getHeader());
		aligned.setReadBases(fqr.getReadBases());
		aligned.setReadName(fqr.getReadName());
		aligned.setBaseQualities(fqr.getBaseQualities());
		aligned.setReferenceIndex(0);
		aligned.setAlignmentStart(0);
		aligned.setCigarString("10S40M50S");
		srfe = new SplitReadFastqExtractor(true, 1, 0, false, false, true, new StringEvidenceIdentifierGenerator());
		List<FastqRecord> result = srfe.extract(aligned);
		Assert.assertEquals(2, result.size());
		Assert.assertEquals(10, result.get(0).getReadLength());
		Assert.assertEquals(50, result.get(1).getReadLength());
	}
}
