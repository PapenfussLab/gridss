package au.edu.wehi.idsv.kraken;

import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.util.Log;
import org.apache.commons.lang.StringUtils;

import java.io.Flushable;
import java.io.IOException;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.concurrent.LinkedBlockingDeque;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Annotates breakend/breakpoint sequences with NCBI taxonomic identifiers.
 * krakenInput and krakenOutput must be pipes connecting to Kraken2 input and output files specified on the command line.
 */
public abstract class KrakenAnnotate<T, U> implements Iterator<U> {
    private static final Log log = Log.getInstance(KrakenAnnotate.class);
    private final OutstandingRecord<T> END_OF_STREAM = new OutstandingRecord<T>(null, null);
    private final Thread feedingThread;
    private final int minSequenceLength;
    private LinkedBlockingDeque<OutstandingRecord<T>> outstanding = new LinkedBlockingDeque<>();
    private KrakenParser krakenOutput;
    private FastqWriter krakenInput;
    private OutstandingRecord<T> nextRecord = null;
    private AtomicInteger recordsInKrakenPipeline = new AtomicInteger(0);

    private static class OutstandingRecord<T> {
        public OutstandingRecord(FastqRecord querySentToKraken, T record) {
            this.querySentToKraken = querySentToKraken;
            this.record = record;
        }
        public final FastqRecord querySentToKraken;
        public final T record;
    }

    /**
     *
     * @param krakenOutput kraken2 output parser.
     *                     The underlying stream is not closed at end of iterator.
     * @param krakenInput Fastq writer connected to kraken2 input.
     *                    Flushed if implements Flushable.
     * @param input VCF records to annotate
     */
    public KrakenAnnotate(
            KrakenParser krakenOutput,
            FastqWriter krakenInput,
            Iterator<T> input,
            int minSequenceLength) {
        this.krakenOutput = krakenOutput;
        this.krakenInput = krakenInput;
        this.minSequenceLength = minSequenceLength;
        this.feedingThread = new Thread(() -> feedKrakenInput(input));
        this.feedingThread.setName("feeding the kraken");
        this.feedingThread.setDaemon(true);
        this.feedingThread.start();
    }

    private void ensureNextRecord() {
        if (nextRecord == null) {
            try {
                nextRecord = outstanding.take();
            } catch (InterruptedException e) {
                log.warn(e);
            }
        }
    }

    protected abstract U transform(T record, KrakenClassification kc);

    @Override
    public boolean hasNext() {
        ensureNextRecord();
        return nextRecord != END_OF_STREAM;
    }

    @Override
    public U next() {
        ensureNextRecord();
        if (nextRecord == END_OF_STREAM) throw new NoSuchElementException();
        U result;
        KrakenClassification kc = null;
        if (nextRecord.querySentToKraken != null) {
            String expectedSequenceId = nextRecord.querySentToKraken.getReadName();
            kc = krakenOutput.next();
            log.debug(String.format("Kraken response: %s %d", kc.sequenceId, kc.taxonomyId));
            if (!expectedSequenceId.equals(kc.sequenceId)) {
                String msg = String.format("Unexpected output from kraken2. Expected result for '%s', found '%s'", expectedSequenceId, kc.sequenceId);
                log.error(msg);
                throw new RuntimeException(msg);
            }
        }
        result = transform(nextRecord.record, kc);
        nextRecord = null;
        return result;
    }

    /**
     * Fastq record to send to Kraken2
     * @param vc variant
     * @return FastqRecord to send to Kraken2, null is no lookup is required.
     */
    protected abstract FastqRecord getKrakenFastqRecord(T vc);

    protected abstract boolean shouldReturnUnprocessedRecords();

    private void process(T record) {
        FastqRecord fqr = getKrakenFastqRecord(record);
        if (fqr != null && fqr.getReadLength() >= minSequenceLength) {
            outstanding.add(new OutstandingRecord(fqr, record));
            recordsInKrakenPipeline.incrementAndGet();
            krakenInput.write(fqr);
        } else {
            if (shouldReturnUnprocessedRecords()) {
                outstanding.add(new OutstandingRecord(null, record));
            }
        }
    }
    private void feedKrakenInput(Iterator<T> input) {
        log.debug("Kraken worker thread started");
        while (input.hasNext()) {
            T record = input.next();
            process(record);
        }
        outstanding.add(END_OF_STREAM);
        krakenInput.close();
        log.debug("Kraken worker thread complete");
    }
}
