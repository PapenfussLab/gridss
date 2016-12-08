package au.edu.wehi.idsv.picard;

import java.util.List;

import org.junit.Assert;
import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;
import htsjdk.samtools.QueryInterval;

public class ReferenceLookupTest extends TestHelper {
	@Test
	public void getIntervals_should_return_valid_contig_intervals() {
		QueryInterval qi = SMALL_FA.getIntervals(1000000).get(0);
		Assert.assertEquals(0, qi.referenceIndex);
		Assert.assertEquals(1, qi.start);
		Assert.assertEquals(10000, qi.end);
	}
	@Test
	public void getIntervals_should_chunk_at_given_size() {
		List<QueryInterval> intervals = SMALL_FA.getIntervals(10);
		Assert.assertEquals(0, intervals.get(0).referenceIndex);
		Assert.assertEquals(1, intervals.get(0).start);
		Assert.assertEquals(10, intervals.get(0).end);
		Assert.assertEquals(0, intervals.get(1).referenceIndex);
		Assert.assertEquals(11, intervals.get(1).start);
		Assert.assertEquals(20, intervals.get(1).end);
		Assert.assertEquals(0, intervals.get(2).referenceIndex);
		Assert.assertEquals(21, intervals.get(2).start);
		Assert.assertEquals(30, intervals.get(2).end);
	}
}
