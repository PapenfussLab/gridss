package gridss.filter;

import au.edu.wehi.idsv.TestHelper;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.SamRecordFilter;
import org.junit.Test;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;


public class UnionAggregateFilterTest extends TestHelper {
	@Test
	public void should_use_underlying_filter() {
		assertTrue(new UnionAggregateFilter(ImmutableList.of(new FilterAll())).filterOut(null));
		assertFalse(new UnionAggregateFilter(ImmutableList.of(new FilterNone())).filterOut(null));
	}
	@Test
	public void should_apply_or_operation() {
		assertFalse(new UnionAggregateFilter(ImmutableList.of(new FilterAll(), new FilterNone())).filterOut(null));
		assertFalse(new UnionAggregateFilter(ImmutableList.of(new FilterNone(), new FilterNone())).filterOut(null));
		assertFalse(new UnionAggregateFilter(ImmutableList.of(new FilterNone(), new FilterAll())).filterOut(null));
		assertTrue(new UnionAggregateFilter(ImmutableList.of(new FilterAll(), new FilterAll())).filterOut(null));
	}
	public static class FilterAll implements SamRecordFilter {
		@Override
		public boolean filterOut(SAMRecord record) {
			return true;
		}
		@Override
		public boolean filterOut(SAMRecord first, SAMRecord second) {
			return true;
		}	
	}
	public static class FilterNone implements SamRecordFilter {
		@Override
		public boolean filterOut(SAMRecord record) {
			return false;
		}
		@Override
		public boolean filterOut(SAMRecord first, SAMRecord second) {
			return false;
		}	
	}
}