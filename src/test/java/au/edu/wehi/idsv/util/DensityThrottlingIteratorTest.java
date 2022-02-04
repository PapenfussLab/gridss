package au.edu.wehi.idsv.util;

import com.google.common.collect.Lists;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import static org.junit.Assert.assertEquals;


public class DensityThrottlingIteratorTest {
	public class IntDensityThrottlingIterator extends DensityThrottlingIterator<Integer> {
		public IntDensityThrottlingIterator(Iterator<Integer> it, int windowSize, double acceptDensity, double targetDensity) {
			super(it, windowSize, acceptDensity, targetDensity);
		}
		@Override
		protected long getPosition(Integer record) {
			return record;
		}

		@Override
		protected boolean excludedFromThrottling(Integer record) {
			return false;
		}
	}
	@Test
	public void should_accept_all_up_to_accept_threshold() {
		List<Integer> input = new ArrayList<Integer>();
		for (int i = 0; i < 1024; i++) {
			input.add(i);
		}
		List<Integer> result = Lists.newArrayList(new IntDensityThrottlingIterator(input.iterator(), 4, 1, 1));
		assertEquals(1024, result.size());
	}
	@Test
	public void should_use_exponential_backoff_to_throttle_above_unconditional_accept_threshold() {
		List<Integer> input = new ArrayList<Integer>();
		for (int i = 0; i < 1024; i++) {
			for (int j = 0; j < 32; j++) {
				input.add(i);
			}
		}
		List<Integer> result = Lists.newArrayList(new IntDensityThrottlingIterator(input.iterator(), 4, 2.0, 4.0));
		assertEquals(2048, result.size(), 64);
	}
	public class IntRecord {
		public final int value;
		public final long position;
		public IntRecord(int value, long position) {
			this.value = value;
			this.position = position;
		}
	}
	public class IntRecordRecordDensityThrottlingIterator extends DensityThrottlingIterator<IntRecord> {
		public IntRecordRecordDensityThrottlingIterator(Iterator<IntRecord> it, int windowSize, double acceptDensity, double targetDensity) {
			super(it, windowSize, acceptDensity, targetDensity);
		}
		@Override
		protected long getPosition(IntRecord record) {
			return record.position;
		}

		@Override
		protected boolean excludedFromThrottling(IntRecord record) {
			return false;
		}
	}
	@Test
	public void should_be_deterministic() {
		List<IntRecord> input = new ArrayList<IntRecord>();
		for (int i = 0; i < 1024; i++) {
			int value = ~i & 0xFF;
			for (long pos = 0; pos < 128; pos++) {
				input.add(new IntRecord(value, pos));
			}
		}
		List<IntRecord> result1 = Lists.newArrayList(new IntRecordRecordDensityThrottlingIterator(input.iterator(), 4, 2.0, 4.0));
		List<IntRecord> result2 = Lists.newArrayList(new IntRecordRecordDensityThrottlingIterator(input.iterator(), 4, 2.0, 4.0));

		assertEquals(result1, result2);
	}
}
