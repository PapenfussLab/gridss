package au.edu.wehi.idsv;

import static org.junit.Assert.*;

import org.junit.Test;

import com.google.common.collect.Iterators;

import htsjdk.samtools.SAMRecord;


public class SpannedIndelEvidenceTest extends TestHelper {	
	@Test
	public void qual_should_use_library_indel_distribution() {
		SAMRecord r = Read(2, 1, "5M5D5M");
		r.setMappingQuality(40);
		SoftClipEvidenceIterator it = new SoftClipEvidenceIterator(SES(), Iterators.singletonIterator(r));
		SpannedIndelEvidence e = (SpannedIndelEvidence)it.next();
		it.close();
		assertEquals(5, e.getBreakendQual(), 0.01);
		assertEquals(5, e.getBreakpointQual(), 0.01);
	}
}
