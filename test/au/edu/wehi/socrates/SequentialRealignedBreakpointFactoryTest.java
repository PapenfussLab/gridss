package au.edu.wehi.socrates;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.List;

import htsjdk.samtools.SAMRecord;

import org.junit.Test;

import com.google.common.collect.Iterators;

public class SequentialRealignedBreakpointFactoryTest extends TestHelper {
	public SequentialRealignedBreakpointFactory getFactory(List<SAMRecord> data) {
		return new SequentialRealignedBreakpointFactory(Iterators.peekingIterator(data.iterator()));
	}
	public class TestDirectedBreakpoint implements DirectedBreakpoint {
		private String id;
		private BreakendSummary location;
		public TestDirectedBreakpoint(int referenceIndex, int start, String id) {
			this.id = id;
			this.location = new BreakendSummary(referenceIndex, BreakendDirection.Forward, start, start, 0);
		}
		@Override
		public BreakendSummary getBreakendSummary() {
			return location;
		}
		@Override
		public String getEvidenceID() {
			return id;
		}
		@Override
		public byte[] getBreakpointSequence() {
			return null;
		}
		@Override
		public byte[] getBreakpointQuality() {
			return null;
		}
	}
	@Test
	public void should_match_by_read_name() {
		SequentialRealignedBreakpointFactory factory = getFactory(toList(withReadName("0#1#n1", Read(0, 1, 1))));
		SAMRecord r = factory.findRealignedSAMRecord(new TestDirectedBreakpoint(0, 1, "n1"));
		assertNotNull(r);
	}
	@Test
	public void should_match_sequential_records() {
		SequentialRealignedBreakpointFactory factory = getFactory(toList(
			withReadName("0#1#n1", Read(0, 1, "1M")),
			withReadName("0#1#n2", Read(0, 1, "1M")),
			withReadName("0#1#n3", Read(0, 1, "1M")),
			withReadName("0#2#n4", Read(0, 1, "1M")),
			withReadName("1#1#n5", Read(0, 1, "1M"))
			));
		assertNotNull(factory.findRealignedSAMRecord(new TestDirectedBreakpoint(0, 1, "n1")));
		assertNotNull(factory.findRealignedSAMRecord(new TestDirectedBreakpoint(0, 1, "n3")));
		assertNotNull(factory.findRealignedSAMRecord(new TestDirectedBreakpoint(0, 1, "n2")));
		assertNotNull(factory.findRealignedSAMRecord(new TestDirectedBreakpoint(0, 2, "n4")));
		assertNotNull(factory.findRealignedSAMRecord(new TestDirectedBreakpoint(1, 1, "n5")));
	}
	@Test(expected=IllegalStateException.class)
	public void should_fail_during_non_sequential_traversal() {
		SequentialRealignedBreakpointFactory factory = getFactory(toList(
				withReadName("0#1#n1", Read(0, 1, 1)),
				withReadName("1#1#n5", Read(0, 1, 1))
				));
		assertNotNull(factory.findRealignedSAMRecord(new TestDirectedBreakpoint(1, 1, "n5")));
		assertNull(factory.findRealignedSAMRecord(new TestDirectedBreakpoint(0, 1, "n1")));
	}
}
