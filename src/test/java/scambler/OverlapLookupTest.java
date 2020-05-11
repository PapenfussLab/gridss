package scambler;

import htsjdk.samtools.SAMRecord;
import org.junit.Test;

import java.util.List;
import java.util.stream.Collectors;

import static org.junit.Assert.assertEquals;

public class OverlapLookupTest extends StringGraphTestHelper {
	@Test
	public void should_report_overlaps() {
		OverlapLookup lookup = new OverlapLookup(5);
		List<SAMRecord> list = overlapping(3, 30, 10);
		List<Read> reads = list.stream()
				.map(r -> Read.create(getContext().getLinear(), r))
				.map(r -> { lookup.add(r); return r; })
				.collect(Collectors.toList());
		assertEquals(2, lookup.successors(reads.get(0)).size());
	}
}
