package au.edu.wehi.idsv.alignment;

import com.google.common.util.concurrent.ThreadFactoryBuilder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.Log;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Runs bwa mem through a JNI interface.
 */
public class BwaStreamingAligner implements StreamingAligner, Closeable {
    private static final Log log = Log.getInstance(BwaStreamingAligner.class);
    private ExecutorService bwaDriver = Executors.newSingleThreadExecutor(new ThreadFactoryBuilder().setDaemon(false).setNameFormat("bwaDriver").build());
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
     * @param bufferSizeInBases number of base pairs of sequence to buffer.
     *                          This buffer is evenly split across the input buffer and buffer to run to bwa.
     *                          Actual invocations to bwa will be with a buffer half this size.
     */
    public BwaStreamingAligner(File reference, SAMSequenceDictionary dict, int threads, int bufferSizeInBases) {
        this.bwaInputBuffer = new LinkedBlockingDeque<>();
        this.aligner = new BwaAligner(reference, dict, threads);
        this.batchSize = bufferSizeInBases / 2 + 1;
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
        if (usedBuffer >= batchSize) {
            processInput();
        }
    }

    // synchronized to ensure record ordering is stable
    private synchronized Future<List<SAMRecord>> processInput() {
        final ArrayList<FastqRecord> inFlightBuffer = new ArrayList<>(bwaInputBuffer.size() + 16);
        int basesSent = 0;
        while (!bwaInputBuffer.isEmpty() && basesSent < batchSize) {
            FastqRecord fq = bwaInputBuffer.poll();
            inFlightBuffer.add(fq);
            basesSent += fq.getReadBases().length;
        }
        if (inFlightBuffer.size() > 0) {
            final int actualBasesSent = basesSent;
            Future<List<SAMRecord>> result = bwaDriver.submit(() -> {
                List<SAMRecord> results = getAligner().align(inFlightBuffer);
                bwaOutputBuffer.addAll(results);
                outstandingBases.addAndGet(-actualBasesSent);
                outstandingRecords.addAndGet(-inFlightBuffer.size());
                return results;
            });
            return result;
        }
        return null;
    }

    @Override
    public void flush() {
        Future<List<SAMRecord>> future = processInput();
        if (future != null) {
            try {
                future.get();
            } catch (InterruptedException e) {
                log.error(e, "Exception flushing bwa results.");
                throw new RuntimeException(e);
            } catch (ExecutionException e) {
                log.error(e, "Exception flushing bwa results.");
                throw new RuntimeException(e);
            }
        }
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
        flush();
        this.bwaDriver.shutdown();
        this.aligner.close();
    }
}
