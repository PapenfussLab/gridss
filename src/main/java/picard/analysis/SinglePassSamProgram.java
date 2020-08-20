/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute, 2019 Daniel Cameron
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.analysis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.argumentcollections.OutputArgumentCollection;
import picard.cmdline.argumentcollections.RequiredOutputArgumentCollection;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.TimeUnit;

/**
 * Super class that is designed to provide some consistent structure between subclasses that
 * simply iterate once over a coordinate sorted BAM and collect information from the records
 * as the go in order to produce some kind of output.
 *
 * @author Tim Fennell
 */
public abstract class SinglePassSamProgram extends CommandLineProgram {
    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input SAM or BAM file.")
    public File INPUT;

    @ArgumentCollection
    public OutputArgumentCollection output = getOutputArgumentCollection();

    protected OutputArgumentCollection getOutputArgumentCollection(){
        return new RequiredOutputArgumentCollection();
    }

    protected File OUTPUT;

    @Argument(doc = "If true (default), then the sort order in the header file will be ignored.",
            shortName = StandardOptionDefinitions.ASSUME_SORTED_SHORT_NAME)
    public boolean ASSUME_SORTED = true;

    @Argument(doc = "Stop after processing N reads, mainly for debugging.")
    public long STOP_AFTER = 0;

    @Argument(doc = "Stop after processing N bases, mainly for debugging.")
    public long STOP_AFTER_BASES = 0;

    private static final Log log = Log.getInstance(SinglePassSamProgram.class);

    @Argument(doc = "Allocate each metrics program it's own thread. I/O and record parsing is still shared.")
    public boolean PROCESS_IN_PARALLEL = true;

    /**
     * Number of SAMRecords to batch together before allocating to worker threads.
     * A larger batch size reduces thread synchronisation overhead.
     */
    private static final int BATCH_SIZE = 512;

    /**
     * Maximum number of outstanding batches.
     * The default of two batches results in double buffering: one batch for I/O and parsing, and another
     * being processed in parallel by the program worker threads.
     */
    private static final int IN_FLIGHT_BATCHES = 2;

    /**
     * End of stream sentinel value.
     * Program worker threads use this object to indicate a unexceptional end of stream.
     */
    private static final Exception EOS_SENTINEL = new Exception();

    /**
     * Set the reference File.
     */
    public void setReferenceSequence(final File referenceFile) {
        REFERENCE_SEQUENCE = referenceFile;
    };

    /**
     * Final implementation of doWork() that checks and loads the input and optionally reference
     * sequence files and the runs the sublcass through the setup() acceptRead() and finish() steps.
     */
    @Override
    protected final int doWork() {
        makeItSo(INPUT, REFERENCE_SEQUENCE, ASSUME_SORTED, STOP_AFTER, STOP_AFTER_BASES, Arrays.asList(this), PROCESS_IN_PARALLEL, PROCESS_IN_PARALLEL);
        return 0;
    }
    public static void makeItSo(final File input,
                                final File referenceSequence,
                                final boolean assumeSorted,
                                final long stopAfter,
                                final Collection<SinglePassSamProgram> programs) {
        makeItSo(input, referenceSequence, assumeSorted, stopAfter, 0, programs);
    }
    public static void makeItSo(final File input,
                                final File referenceSequence,
                                final boolean assumeSorted,
                                final long stopAfter,
                                final long stopAfterBases,
                                final Collection<SinglePassSamProgram> programs) {
        makeItSo(input, referenceSequence, assumeSorted, stopAfter, stopAfterBases, programs, true, true);
    }
    public static void makeItSo(final File input,
                                final File referenceSequence,
                                final boolean assumeSorted,
                                final long stopAfter,
                                final long stopAfterBases,
                                final Collection<SinglePassSamProgram> programs,
                                boolean parallel,
                                boolean useAsyncIterator) {

        // Setup the standard inputs
        IOUtil.assertFileIsReadable(input);
        final SamReader in = SamReaderFactory.makeDefault().referenceSequence(referenceSequence).open(input);

        // Optionally load up the reference sequence and double check sequence dictionaries
        final ReferenceSequenceFileWalker walker;
        if (referenceSequence == null) {
            walker = null;
        } else {
            IOUtil.assertFileIsReadable(referenceSequence);
            walker = new ReferenceSequenceFileWalker(referenceSequence);

            if (!in.getFileHeader().getSequenceDictionary().isEmpty()) {
                SequenceUtil.assertSequenceDictionariesEqual(in.getFileHeader().getSequenceDictionary(),
                        walker.getSequenceDictionary());
            }
        }

        // Check on the sort order of the BAM file
        {
            final SortOrder sort = in.getFileHeader().getSortOrder();
            if (sort != SortOrder.coordinate) {
                if (assumeSorted) {
                    log.warn("File reports sort order '" + sort + "', assuming it's coordinate sorted anyway.");
                } else {
                    throw new PicardException("File " + input.getAbsolutePath() + " should be coordinate sorted but " +
                            "the header says the sort order is " + sort + ". If you believe the file " +
                            "to be coordinate sorted you may pass ASSUME_SORTED=true");
                }
            }
        }

        final List<ArrayBlockingQueue<List<Tuple<ReferenceSequence, SAMRecord>>>> buffers = new ArrayList<>(programs.size());
        final List<SinglePassSamProgramRunner> workers = new ArrayList<>(programs.size());
        // Call the abstract setup method!
        boolean anyUseNoRefReads = false;
        for (final SinglePassSamProgram program : programs) {
            if (program.OUTPUT == null) {
                program.OUTPUT = program.output.getOutputFile();
            }
            program.setup(in.getFileHeader(), input);
            anyUseNoRefReads = anyUseNoRefReads || program.usesNoRefReads();

            if (parallel) {
                ArrayBlockingQueue<List<Tuple<ReferenceSequence, SAMRecord>>> buffer = new ArrayBlockingQueue<>(IN_FLIGHT_BATCHES);
                SinglePassSamProgramRunner runner = new SinglePassSamProgramRunner(program, buffer);
                buffers.add(buffer);
                workers.add(runner);
                Thread t = new Thread(runner);
                t.setName(program.toString());
                t.setDaemon(true);
                t.start();
            }
        }

        final ProgressLogger progress = new ProgressLogger(log, 10000000);
        try (CloseableIterator<SAMRecord> it = useAsyncIterator ? new AsyncBufferedIterator<>(in.iterator(), BATCH_SIZE, IN_FLIGHT_BATCHES, "SinglePassSamProgram") : in.iterator()){
            int basesProcessed = 0;
            List<Tuple<ReferenceSequence, SAMRecord>> batch = new ArrayList<>();
            while (it.hasNext()) {
                final SAMRecord rec =  it.next();
                final ReferenceSequence ref;
                if (walker == null || rec.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
                    ref = null;
                } else {
                    ref = walker.get(rec.getReferenceIndex());
                }

                if (parallel) {
                    batch.add(new Tuple<>(ref, rec));
                    if (batch.size() >= BATCH_SIZE) {
                        asyncAcceptReads(buffers, workers, batch);
                        batch = new ArrayList<>();
                    }
                } else {
                    for (final SinglePassSamProgram program : programs) {
                        program.acceptRead(rec, ref);
                    }
                }

                progress.record(rec);
                basesProcessed += rec.getReadLength();

                // See if we need to terminate early?
                if (stopAfter > 0 && progress.getCount() >= stopAfter) {
                    break;
                }
                if (stopAfterBases > 0 && basesProcessed >= stopAfterBases) {
                    break;
                }

                // And see if we're into the unmapped reads at the end
                if (!anyUseNoRefReads && rec.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
                    break;
                }
            }
            if (parallel) {
                if (batch.size() > 0) {
                    asyncAcceptReads(buffers, workers, batch);
                }
                asyncAcceptReads(buffers, workers, new ArrayList<>()); // Empty batch is the EOS indicator
                asyncWaitForCompletion(workers);
            } else {
                for (final SinglePassSamProgram program : programs) {
                    program.finish();
                }
            }
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        } finally {
            CloserUtil.close(in);
        }
    }
    private static void asyncAcceptReads(
            final List<ArrayBlockingQueue<List<Tuple<ReferenceSequence, SAMRecord>>>> buffers,
            final List<SinglePassSamProgramRunner> workers,
            final List<Tuple<ReferenceSequence, SAMRecord>> batch) throws InterruptedException {
        for (int i = 0; i < workers.size(); i++) {
            asyncAcceptRead(buffers.get(i), workers.get(i), batch);
        }
    }
    private static void asyncAcceptRead(
            final ArrayBlockingQueue<List<Tuple<ReferenceSequence, SAMRecord>>> buffer,
            final SinglePassSamProgramRunner worker,
            final List<Tuple<ReferenceSequence, SAMRecord>> batch) throws InterruptedException {
        // Propagate exceptions on worker threads back to main
        while (!buffer.offer(batch, 1, TimeUnit.SECONDS)) {
            // Check if the worker thread is still alive
            raiseAsyncException(worker);
            if (worker.isComplete()) {
                throw new RuntimeException(worker.program.getClass().getName() + " terminated before all records read.");
            }
        }
    }
    private static void raiseAsyncException(final SinglePassSamProgramRunner worker) {
        Exception e = worker.getException();
        if (e != null) {
            throw new RuntimeException("Exception when running " + worker.program.getClass().getName(), e);
        }
    }
    private static void asyncWaitForCompletion(final List<SinglePassSamProgramRunner> workers) throws InterruptedException {
        for (SinglePassSamProgramRunner worker : workers) {
            asyncWaitForCompletion(worker);
        }
    }
    private static void asyncWaitForCompletion(final SinglePassSamProgramRunner worker) throws InterruptedException {
        while (!worker.isComplete()) {
            raiseAsyncException(worker);
            Thread.sleep(50);
        }
    }
    private static class SinglePassSamProgramRunner implements Runnable {
        private final ArrayBlockingQueue<List<Tuple<ReferenceSequence, SAMRecord>>> buffer;
        private final SinglePassSamProgram program;
        private volatile boolean isComplete = false;
        private volatile Exception exception = null;
        public SinglePassSamProgramRunner(SinglePassSamProgram program,  ArrayBlockingQueue<List<Tuple<ReferenceSequence, SAMRecord>>> buffer) {
            this.program = program;
            this.buffer = buffer;
        }

        public boolean isComplete() {
            return isComplete;
        }

        public Exception getException() {
            return exception;
        }

        @Override
        public void run() {
            try {
                List<Tuple<ReferenceSequence, SAMRecord>> batch = buffer.take();
                while (batch.size() > 0) {
                    for (Tuple<ReferenceSequence, SAMRecord> r : batch) {
                        program.acceptRead(r.b, r.a);
                    }
                    batch = buffer.take();
                }
                program.finish();
            } catch (Exception e) {
                exception = e;
            } finally {
                isComplete = true;
            }
        }
    }

    /** Can be overridden and set to false if the section of unmapped reads at the end of the file isn't needed. */
    protected boolean usesNoRefReads() { return true; }

    /** Should be implemented by subclasses to do one-time initialization work. */
    protected abstract void setup(final SAMFileHeader header, final File samFile);

    /**
     * Should be implemented by subclasses to accept SAMRecords one at a time.
     * If the read has a reference sequence and a reference sequence file was supplied to the program
     * it will be passed as 'ref'. Otherwise 'ref' may be null.
     */
    protected abstract void acceptRead(final SAMRecord rec, final ReferenceSequence ref);

    /** Should be implemented by subclasses to do one-time finalization work. */
    protected abstract void finish();

}
