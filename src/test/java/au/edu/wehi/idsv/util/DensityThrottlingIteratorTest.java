package au.edu.wehi.idsv.util;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.Lists;


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
}
