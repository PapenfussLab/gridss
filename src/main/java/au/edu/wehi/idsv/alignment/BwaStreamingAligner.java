package au.edu.wehi.idsv.alignment;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.LinkedBlockingDeque;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Runs bwa mem through a JNI interface.
 */
public class BwaStreamingAligner implements StreamingAligner {
    private static final Log log = Log.getInstance(BwaStreamingAligner.class);
    private Queue<FastqRecord> bwaInputBuffer;
    private final Queue<SAMRecord> bwaOutputBuffer = new LinkedBlockingDeque<>();
    private final BwaAligner aligner;
    private AtomicInteger outstanding = new AtomicInteger(0);
    private AtomicInteger processed = new AtomicInteger(0);
    private AtomicInteger total = new AtomicInteger(0);

    public BwaAligner getAligner() {
        return this.aligner;
    }

    public BwaStreamingAligner(File reference, SAMSequenceDictionary dict, int threads, int bufferSize) {
        this.bwaInputBuffer = new ArrayBlockingQueue<>(bufferSize);
        this.aligner = new BwaAligner(reference, dict, threads);
    }

    @Override
    public void asyncAlign(FastqRecord fq) {
        bwaInputBuffer.add(fq);
        outstanding.incrementAndGet();
        total.incrementAndGet();
    }

    private synchronized void processInput() {
        ArrayList<FastqRecord> inFlightBuffer = new ArrayList<>();
        while (!bwaInputBuffer.isEmpty()) {
            // Technically we have a race condition here as we could add records faster
            // than we can buffer them thus causing our buffer size to be unbounded
            inFlightBuffer.add(bwaInputBuffer.poll());
        }
        if (inFlightBuffer.size() > 0) {
            List<SAMRecord> results = aligner.align(inFlightBuffer);
            bwaOutputBuffer.addAll(results);
        }
    }

    @Override
    public void flush() throws IOException {
        processInput();
    }

    @Override
    public boolean hasAlignmentRecord() {
        return total.get() > 0;
    }

    @Override
    public int processedAlignmentRecords() {
        return processed.get();
    }

    @Override
    public int outstandingAlignmentRecord() {
        return outstanding.get();
    }

    @Override
    public SAMRecord getAlignment() {
        if (total.decrementAndGet() < 0) {
            total.incrementAndGet(); // We're not actually consuming a record - put it back
            throw new IllegalStateException("No more records left to process");
        }
        if (bwaOutputBuffer.isEmpty()) {
            processInput();
        }
        return bwaOutputBuffer.poll();
    }

    @Override
    public void close() throws IOException {
        this.aligner.close();
    }
}
