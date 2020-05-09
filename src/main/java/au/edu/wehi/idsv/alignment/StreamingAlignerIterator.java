package au.edu.wehi.idsv.alignment;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.Log;

import java.io.IOException;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.concurrent.atomic.AtomicBoolean;

public class StreamingAlignerIterator implements StreamingAligner, Iterator<SAMRecord> {
    private static final Log log = Log.getInstance(StreamingAlignerIterator.class);
    private static final int POLL_INTERVAL = 100;
    private final StreamingAligner underlying;
    private final AtomicBoolean isClosed = new AtomicBoolean(false);

    public StreamingAlignerIterator(StreamingAligner underlying) {
        this.underlying = underlying;
    }
    @Override
    public void asyncAlign(FastqRecord fq) throws IOException {
        underlying.asyncAlign(fq);
    }

    /**
     * Flushes outstanding alignments.
     */
    @Override
    public void flush() throws IOException {
        underlying.flush();
    }

    @Override
    public int processedAlignmentRecords() {
        return underlying.processedAlignmentRecords();
    }

    @Override
    public int outstandingAlignmentRecord() {
        return underlying.outstandingAlignmentRecord();
    }

    @Override
    public SAMRecord getAlignment() {
        return underlying.getAlignment();
    }

    /**
     * Flushes outstanding alignments.
     * Alignment records returned by the aligner are still available after closing.
     */
    @Override
    public synchronized void close() throws IOException {
        flush();
        isClosed.set(true);
    }

    private void syncEnsureNext() {
        while (!isClosed.get() && processedAlignmentRecords() == 0) {
            try {
                log.debug(String.format("%d alignments outstanding", outstandingAlignmentRecord()));
                Thread.sleep(POLL_INTERVAL);
            } catch (InterruptedException e) {
                log.warn(e);
                return;
            }
        }
    }
    /**
     * Blocks until it the next alignment record is available.
     */
    @Override
    public boolean hasNext() {
        syncEnsureNext();
        return processedAlignmentRecords() > 0;
    }
    /**
     * Blocks until the aligner returns the next alignment.
     */
    @Override
    public SAMRecord next() {
        if (!hasNext()) throw new NoSuchElementException();
        return getAlignment();
    }
}
