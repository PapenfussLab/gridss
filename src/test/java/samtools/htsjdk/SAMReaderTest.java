package samtools.htsjdk;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReader.AssertingIterator;
import htsjdk.samtools.util.DelegatingIterator;
import org.junit.Test;

import java.util.ArrayList;

public class SAMReaderTest {
	@Test(expected=IllegalStateException.class)
	//public void AssertingIterator_should_not_retain_reference_to_emitted_records() {
	public void when_this_starts_failing_remove_workaround_SAMEvidenceSource_transform_workaround() {
		SAMFileHeader header = new SAMFileHeader();
		header.addSequence(new SAMSequenceRecord("contig1", 10));
		header.addSequence(new SAMSequenceRecord("contig2", 10));
		SAMRecord r1 = new SAMRecord(header);
		r1.setReferenceIndex(0);
		r1.setAlignmentStart(1);
		SAMRecord r2 = new SAMRecord(header);
		r2.setReferenceIndex(0);
		r2.setAlignmentStart(2);
		ArrayList<SAMRecord> list = new ArrayList<>();
		list.add(r1);
		list.add(r2);
		AssertingIterator it = new SamReader.AssertingIterator(new DelegatingIterator<SAMRecord>(list.iterator()));
		it.assertSorted(SortOrder.coordinate);
		it.next();
		// the iterator has already returned r1. It shouldn't matter if I now change it
		// it was in the correct order when it got iterated over
		r1.setAlignmentStart(3);
		r1.setReferenceIndex(1);
		r1.setReadUnmappedFlag(true);
		it.next();
		it.close();
	}
}
