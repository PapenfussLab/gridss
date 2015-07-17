package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.ProgressLoggerInterface;

import java.util.Iterator;

public class ProgressLoggingSAMRecordIterator implements CloseableIterator<SAMRecord> {
	private final ProgressLoggerInterface logger;
	private final Iterator<SAMRecord> iterator;
	public ProgressLoggingSAMRecordIterator(Iterator<SAMRecord> iterator, ProgressLoggerInterface logger) {
		this.iterator = iterator;
		this.logger = logger;
	}
	@Override
	public boolean hasNext() {
		return iterator.hasNext();
	}
	@Override
	public SAMRecord next() {
		SAMRecord n = iterator.next();
		logger.record(n);
		return n;
	}
	@Override
	public void close() {
		CloserUtil.close(iterator);
	}
}
