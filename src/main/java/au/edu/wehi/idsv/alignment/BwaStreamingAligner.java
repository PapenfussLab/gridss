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
    private final int batchSize;
    private Queue<FastqRecord> bwaInputBuffer;
    private final Queue<SAMRecord> bwaOutputBuffer = new LinkedBlockingDeque<>();
    private final BwaAligner aligner;
    private AtomicInteger outstandingRecords = new AtomicInteger(0);
    private AtomicInteger outstandingBases = new AtomicInteger(0);
    public BwaAligner getAligner() {
        return this.aligner;
    }

    /**
     *
     * @param reference Reference genome
     * @param dict sequence dictionary for reference genome
     * @param threads number of bwa threads
     * @param batchSize number of base pairs of sequence to batch before calling BWA
     */
    public BwaStreamingAligner(File reference, SAMSequenceDictionary dict, int threads, int batchSize) {
        this.bwaInputBuffer = new LinkedBlockingDeque<>();
        this.aligner = new BwaAligner(reference, dict, threads);
        this.batchSize = batchSize;
    }

    /**
     * Align the given records.
     *
     * @implNote This method currently uses a buffered synchronous implementation.
     * @param fq
     */
    @Override
    public void asyncAlign(FastqRecord fq) {
        bwaInputBuffer.add(fq);
        outstandingRecords.incrementAndGet();
        int usedBuffer = outstandingBases.addAndGet(fq.getReadBases().length);
        // TODO: actually make this asynchronous
        if (usedBuffer >= batchSize) {
            processInput();
        }
    }

    private synchronized void processInput() {
        ArrayList<FastqRecord> inFlightBuffer = new ArrayList<>(bwaInputBuffer.size() + 16);
        int basesSent = 0;
        while (!bwaInputBuffer.isEmpty() && basesSent < batchSize) {
            FastqRecord fq = bwaInputBuffer.poll();
            inFlightBuffer.add(fq);
            basesSent += fq.getReadBases().length;
        }
        if (inFlightBuffer.size() > 0) {
            List<SAMRecord> results = aligner.align(inFlightBuffer);
            bwaOutputBuffer.addAll(results);
        }
        outstandingBases.addAndGet(-basesSent);
        outstandingRecords.addAndGet(-inFlightBuffer.size());
    }

    @Override
    public void flush() throws IOException {
        processInput();
    }

    @Override
    public int processedAlignmentRecords() {
        return bwaOutputBuffer.size();
    }

    @Override
    public int outstandingAlignmentRecord() {
        return outstandingRecords.get();
    }

    @Override
    public SAMRecord getAlignment() {
        SAMRecord result = bwaOutputBuffer.poll();
        if (result == null) {
            throw new IllegalStateException("Call flush() or check processedAlignmentRecords() to ensure records are available.");
        }
        return result;
    }

    @Override
    public void close() throws IOException {
        this.aligner.close();
    }
}
