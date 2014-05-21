package au.edu.wehi.socrates;

import static org.junit.Assert.assertEquals;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.ProgressLoggerInterface;

import java.util.Iterator;
import java.util.List;

import org.junit.Test;

import com.google.common.collect.Lists;

public class ProgressLoggingSAMRecordIteratorTest extends TestHelper {
	public class StubProgressLoggerInterface implements ProgressLoggerInterface {
		public int count = 0;
		@Override
		public boolean record(String chrom, int pos) {
			count++;
			return true;
		}
		@Override
		public boolean record(SAMRecord rec) {
			count++;
			return true;
		}
		@Override
		public boolean record(SAMRecord... recs) {
			count += recs.length;
			return true;
		}
	}
	@Test
	public void should_call_logger() {
		StubProgressLoggerInterface log = new StubProgressLoggerInterface();
		List<SAMRecord> list = Lists.newArrayList(RP(1, 2, 3));
		for (Iterator<SAMRecord> it = new ProgressLoggingSAMRecordIterator(list.iterator(), log); it.hasNext(); it.next()) {
		}
		assertEquals(2, log.count);
	}
	@Test
	public void should_allow_null_logger() {
		List<SAMRecord> list = Lists.newArrayList(RP(1, 2, 3));
		for (Iterator<SAMRecord> it = new ProgressLoggingSAMRecordIterator(list.iterator(), null); it.hasNext(); it.next()) {	
		}
	}
}
