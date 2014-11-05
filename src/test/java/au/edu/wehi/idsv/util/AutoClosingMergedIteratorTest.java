package au.edu.wehi.idsv.util;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;

import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;

import com.google.common.collect.ImmutableList;


public class AutoClosingMergedIteratorTest extends TestHelper {
	@Test
	public void should_terminate_async_reader_thread_for_underlying_stream() throws InterruptedException {
		List<SAMRecord> list = ImmutableList.of(Read(0, 1, 1), Read(0, 1, 1), Read(0, 1, 1), Read(0, 1, 1), Read(0, 2, 1));
		AsyncBufferedIterator<SAMRecord> abi = new AsyncBufferedIterator<SAMRecord>(list.iterator(), 1, 1);
		AutoClosingMergedIterator<SAMRecord> merged = new AutoClosingMergedIterator<SAMRecord>(ImmutableList.of(
				abi), new SAMRecordCoordinateComparator());
		merged.next();
		Thread.sleep(1);
		assertNotNull(AsyncBufferedIteratorTest.getThreadWithName(abi.getBackgroundThreadName()));
		merged.close();
		Thread.sleep(1);
		assertNull(AsyncBufferedIteratorTest.getThreadWithName(abi.getBackgroundThreadName()));
	}
}
