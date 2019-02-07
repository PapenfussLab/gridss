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
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import javax.management.InstanceNotFoundException;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Executor;
import java.util.concurrent.LinkedBlockingDeque;
import java.util.concurrent.atomic.AtomicInteger;

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

    @Argument(shortName = "O", doc = "File to write the output to.")
    public File OUTPUT;

    @Argument(doc = "If true (default), then the sort order in the header file will be ignored.",
            shortName = StandardOptionDefinitions.ASSUME_SORTED_SHORT_NAME)
    public boolean ASSUME_SORTED = true;

    @Argument(doc = "Stop after processing N reads, mainly for debugging.")
    public long STOP_AFTER = 0;

    private static final Log log = Log.getInstance(SinglePassSamProgram.class);

    private static final int BATCH_SIZE = 512;
    private static final int IN_FLIGHT_BATCHES = 2;
    private static final Exception EOS_SENTINEL = new Exception();
    private static final boolean USE_ASYNC_ITERATOR = true;
    private static final boolean USE_ASYNC = true;

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
        makeItSo(INPUT, REFERENCE_SEQUENCE, ASSUME_SORTED, STOP_AFTER, Arrays.asList(this));
        return 0;
    }

    public static void makeItSo(final File input,
                                final File referenceSequence,
                                final boolean assumeSorted,
                                final long stopAfter,
                                final Collection<SinglePassSamProgram> programs) {

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

        final BlockingQueue<Exception> completion = new LinkedBlockingDeque<>();
        final List<ArrayBlockingQueue<List<Tuple<ReferenceSequence, SAMRecord>>>> buffers = new ArrayList<>(programs.size());
        // Call the abstract setup method!
        boolean anyUseNoRefReads = false;
        for (final SinglePassSamProgram program : programs) {
            program.setup(in.getFileHeader(), input);
            anyUseNoRefReads = anyUseNoRefReads || program.usesNoRefReads();

            if (USE_ASYNC) {
                ArrayBlockingQueue<List<Tuple<ReferenceSequence, SAMRecord>>> buffer = new ArrayBlockingQueue<>(IN_FLIGHT_BATCHES);
                buffers.add(buffer);
                Thread t = new Thread(new SinglePassSamProgramRunner(program, buffer, completion));
                t.setName(program.toString());
                t.setDaemon(true);
                t.start();
            }
        }

        final ProgressLogger progress = new ProgressLogger(log, 10000000);
        try (CloseableIterator<SAMRecord> it = USE_ASYNC_ITERATOR ? new AsyncBufferedIterator<>(in.iterator(), BATCH_SIZE, IN_FLIGHT_BATCHES, "SinglePassSamProgram") : in.iterator()){
            List<Tuple<ReferenceSequence, SAMRecord>> batch = new ArrayList<>();
            while (it.hasNext()) {
                final SAMRecord rec =  it.next();
                final ReferenceSequence ref;
                if (walker == null || rec.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
                    ref = null;
                } else {
                    ref = walker.get(rec.getReferenceIndex());
                }

                if (USE_ASYNC) {
                    batch.add(new Tuple<>(ref, rec));
                    if (batch.size() >= BATCH_SIZE) {
                        for (ArrayBlockingQueue<List<Tuple<ReferenceSequence, SAMRecord>>> buffer : buffers) {
                            buffer.put(batch);
                        }
                        batch = new ArrayList<>();
                    }
                    // Raise any encountered exceptions as early as possible
                    Exception ex = completion.peek();
                    if (ex != null && ex != EOS_SENTINEL) {
                        if (ex instanceof RuntimeException) {
                            throw (RuntimeException) ex;
                        } else {
                            throw new RuntimeException(ex);
                        }
                    }
                } else {
                    for (final SinglePassSamProgram program : programs) {
                        program.acceptRead(rec, ref);
                    }
                }

                progress.record(rec);

                // See if we need to terminate early?
                if (stopAfter > 0 && progress.getCount() >= stopAfter) {
                    break;
                }

                // And see if we're into the unmapped reads at the end
                if (!anyUseNoRefReads && rec.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
                    break;
                }
            }
            if (USE_ASYNC) {
                batch.add(new Tuple<>(null, null)); // Add EOS indicator
                for (ArrayBlockingQueue<List<Tuple<ReferenceSequence, SAMRecord>>> buffer : buffers) {
                    buffer.put(batch);
                }
                for (int i = 0; i < programs.size(); i++) {
                    Exception ex = completion.take();
                    if (ex != EOS_SENTINEL) {
                        if (ex instanceof RuntimeException) {
                            throw (RuntimeException) ex;
                        } else {
                            throw new RuntimeException(ex);
                        }
                    }
                }
            } else {
                for (final SinglePassSamProgram program : programs) {
                    program.finish();
                }
            }
        } catch (InterruptedException ex) {
            throw new RuntimeException(ex);
        } finally {
            CloserUtil.close(in);
        }
    }
    private static class SinglePassSamProgramRunner implements Runnable {
        private final ArrayBlockingQueue<List<Tuple<ReferenceSequence, SAMRecord>>> buffer;
        private final SinglePassSamProgram program;
        private final BlockingQueue<Exception> completion;
        public SinglePassSamProgramRunner(SinglePassSamProgram program,
                                          ArrayBlockingQueue<List<Tuple<ReferenceSequence, SAMRecord>>> buffer,
                                          BlockingQueue<Exception> completion) {
            this.program = program;
            this.buffer = buffer;
            this.completion = completion;
        }

        @Override
        public void run() {
            Exception exception = EOS_SENTINEL;
            try {
                List<Tuple<ReferenceSequence, SAMRecord>> batch = buffer.take();
                while (true) {
                    for (Tuple<ReferenceSequence, SAMRecord> r : batch) {
                        if (r.b == null) {
                            program.finish();
                            throw EOS_SENTINEL;
                        }
                        program.acceptRead(r.b, r.a);
                    }
                    batch = buffer.take();
                }
            } catch (Exception e) {
                exception = e;
            } finally {
                completion.add(exception);
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