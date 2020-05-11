package au.edu.wehi.idsv.picard;

import au.edu.wehi.idsv.TestHelper;
import htsjdk.samtools.QueryInterval;
import org.junit.Assert;
import org.junit.Test;

import java.util.List;

public class ReferenceLookupTest extends TestHelper {
	@Test
	public void getIntervals_should_return_valid_contig_intervals() {
		QueryInterval[] qi = SMALL_FA.getIntervals(1000000, 1000000).get(0);
		Assert.assertEquals(0, qi[0].referenceIndex);
		Assert.assertEquals(1, qi[0].start);
		Assert.assertEquals(10000, qi[0].end);
	}
	@Test
	public void getIntervals_should_chunk_at_given_size() {
		List<QueryInterval[]> intervals = SMALL_FA.getIntervals(10, 10);
		Assert.assertEquals(0, intervals.get(0)[0].referenceIndex);
		Assert.assertEquals(1, intervals.get(0)[0].start);
		Assert.assertEquals(10, intervals.get(0)[0].end);
		Assert.assertEquals(0, intervals.get(1)[0].referenceIndex);
		Assert.assertEquals(11, intervals.get(1)[0].start);
		Assert.assertEquals(20, intervals.get(1)[0].end);
		Assert.assertEquals(0, intervals.get(2)[0].referenceIndex);
		Assert.assertEquals(21, intervals.get(2)[0].start);
		Assert.assertEquals(30, intervals.get(2)[0].end);
	}
	@Test
	public void getIntervals_should_add_change_penalty_when_moving_to_next_contig() {
		InMemoryReferenceSequenceFile ref = new InMemoryReferenceSequenceFile(new String[] {"chr1", "chr2", "chr3" }, new byte[][] { B("AAAAA"), B("TTTT"), B("GGGG") } );
		// 12345 1234 1234
		// *******       *
		//        *******
		List<QueryInterval[]> intervals = ref.getIntervals(7, 1);
		ref.close();
		Assert.assertEquals(2, intervals.get(0).length);
		Assert.assertEquals(0, intervals.get(0)[0].referenceIndex);
		Assert.assertEquals(1, intervals.get(0)[0].start);
		Assert.assertEquals(5, intervals.get(0)[0].end);
		Assert.assertEquals(1, intervals.get(0)[1].referenceIndex);
		Assert.assertEquals(1, intervals.get(0)[1].start);
		Assert.assertEquals(1, intervals.get(0)[1].end);
		
		Assert.assertEquals(1, intervals.get(1)[0].referenceIndex);
		Assert.assertEquals(2, intervals.get(1)[0].start);
		Assert.assertEquals(4, intervals.get(1)[0].end);
		Assert.assertEquals(2, intervals.get(1)[1].referenceIndex);
		Assert.assertEquals(1, intervals.get(1)[1].start);
		Assert.assertEquals(3, intervals.get(1)[1].end);
		
		Assert.assertEquals(2, intervals.get(2)[0].referenceIndex);
		Assert.assertEquals(4, intervals.get(2)[0].start);
		Assert.assertEquals(4, intervals.get(2)[0].end);
	}
}
