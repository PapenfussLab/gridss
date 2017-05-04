package au.edu.wehi.idsv.util;

import org.junit.Assert;
import org.junit.Test;

import htsjdk.samtools.util.Log;

public class MessageThrottlerTest {
	private static final Log log = Log.getInstance(MessageThrottlerTest.class);
	@Test
	public void should_suppress_after_n_messages() {
		MessageThrottler mt = new MessageThrottler(3);
		Assert.assertFalse(mt.shouldSupress(log, "test message"));
		Assert.assertFalse(mt.shouldSupress(log, "test message"));
		Assert.assertFalse(mt.shouldSupress(log, "test message"));
		Assert.assertTrue(mt.shouldSupress(log, "test message"));
		Assert.assertTrue(mt.shouldSupress(log, "test message"));
	}
}
