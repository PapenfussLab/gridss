package au.edu.wehi.idsv;

import static org.junit.Assert.*;

import org.junit.Test;

import au.edu.wehi.idsv.bed.IntervalBed;


public class SequentialCoverageThresholdTest extends TestHelper {
	@Test
	public void should_calculate_aligned_coverage() {
		SequentialCoverageThreshold sct = new SequentialCoverageThreshold(getSequenceDictionary(), getContext().getLinear(), 2);
		sct.acceptRecord(Read(0, 1, "10M")); // 1-10
		sct.acceptRecord(Read(0, 2, "1S10M5S")); // 2-11
		sct.acceptRecord(Read(1, 1, "3M"));
		sct.acceptRecord(Read(1, 1, "3M"));
		sct.acceptRecord(Read(1, 1, "3M"));
		IntervalBed bed = sct.finish();
		for (int i = 0; i < 100; i++) {
			assertEquals(i >= 2 && i <= 10, bed.overlaps(0, i, i));
			assertEquals(i >= 1 && i <= 3, bed.overlaps(1, i, i));
		}
	}
}
