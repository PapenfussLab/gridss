package au.edu.wehi.validation;

import org.junit.Test;

import java.io.File;

import static org.junit.Assert.assertTrue;


public class LongReadSupportFinderTest {
	@Test
	public void pacbio_spanning_reads_should_sum_deletion_events() {
		File file = new File("src/test/resources/pacbiona12989chem1chr1_196132675-196183463.bam");
		LongReadSupportFinder lr = new LongReadSupportFinder(file);
		LongReadSupportLevel support = lr.evaluateDeletion("chr1", 196158052, 196158052, 196158331, 196158331,
				new ValidateDeletions().SOFT_CLIP_MARGIN, new ValidateDeletions().MIN_SOFT_CLIP_LENGTH, new ValidateDeletions().SPANNING_WINDOW_SIZE, new ValidateDeletions().MIN_DELETION_LENGTH);
		//read m130208_134627_42137_c100474062550000001823070606131326_s1_p0/79081/0_10472
		assertTrue(support.spanningAlignments.contains(272));
	}
}
