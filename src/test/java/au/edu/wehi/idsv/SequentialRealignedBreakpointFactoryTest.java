package au.edu.wehi.idsv;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import htsjdk.samtools.SAMRecord;

import java.util.List;

import org.junit.Test;

import com.google.common.collect.Iterators;

public class SequentialRealignedBreakpointFactoryTest extends TestHelper {
	public SequentialRealignedBreakpointFactory getFactory(List<SAMRecord> data) {
		return new SequentialRealignedBreakpointFactory(Iterators.peekingIterator(data.iterator()));
	}
	@Test
	public void should_match_by_read_name() {
		SequentialRealignedBreakpointFactory factory = getFactory(L(withReadName("0#1#n1", Read(0, 1, 1))));
		SAMRecord r = factory.findAssociatedSAMRecord(new MockDirectedEvidence(0, 1, "n1"));
		assertNotNull(r);
	}
	@Test
	public void should_match_sequential_records() {
		SequentialRealignedBreakpointFactory factory = getFactory(L(
			withReadName("0#1#n1", Read(0, 1, "1M")),
			withReadName("0#1#n2", Read(0, 1, "1M")),
			withReadName("0#1#n3", Read(0, 1, "1M")),
			withReadName("0#2#n4", Read(0, 1, "1M")),
			withReadName("1#1#n5", Read(0, 1, "1M"))
			));
		assertNotNull(factory.findAssociatedSAMRecord(new MockDirectedEvidence(0, 1, "n1")));
		assertNotNull(factory.findAssociatedSAMRecord(new MockDirectedEvidence(0, 1, "n3")));
		assertNotNull(factory.findAssociatedSAMRecord(new MockDirectedEvidence(0, 1, "n2")));
		assertNotNull(factory.findAssociatedSAMRecord(new MockDirectedEvidence(0, 2, "n4")));
		assertNotNull(factory.findAssociatedSAMRecord(new MockDirectedEvidence(1, 1, "n5")));
	}
	@Test(expected=IllegalStateException.class)
	public void should_fail_during_non_sequential_traversal() {
		SequentialRealignedBreakpointFactory factory = getFactory(L(
				withReadName("0#1#n1", Read(0, 1, 1)),
				withReadName("1#1#n5", Read(0, 1, 1))
				));
		assertNotNull(factory.findAssociatedSAMRecord(new MockDirectedEvidence(1, 1, "n5")));
		assertNull(factory.findAssociatedSAMRecord(new MockDirectedEvidence(0, 1, "n1")));
	}
	@Test(expected=RuntimeException.class)
	public void should_error_if_expected_record_not_found() {
		SequentialRealignedBreakpointFactory factory = getFactory(L(
			withReadName("0#1#n1", Read(0, 1, "1M"))
			));
		factory.findAssociatedSAMRecord(new MockDirectedEvidence(0, 1, "n2"), true);
	}
}
