package gridss.cmdline;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Locale;

import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.picard.TwoBitBufferedReferenceSequenceFile;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SequenceUtil;
import picard.PicardException;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

/**
 * Base class used to transform a VCF breakpoint call set given the full evidence available.
 * 
 * 
 * @author Daniel Cameron
 *
 */
public abstract class ByReadNameSinglePassSamProgram extends ReferenceCommandLineProgram {
	private static final Log log = Log.getInstance(ByReadNameSinglePassSamProgram.class);
	@Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input SAM or BAM file. " +
			" If multiple mapping locations are reported for each read, these reads must be grouped together.")
    public File INPUT;

    @Option(shortName = "O", doc = "File to write the output to.")
    public File OUTPUT;

    @Option(doc = "If true (default), then the sort order in the header file will be ignored.",
            shortName = StandardOptionDefinitions.ASSUME_SORTED_SHORT_NAME)
    public boolean ASSUME_SORTED = true;

    @Option(doc = "Stop after processing N reads, mainly for debugging.")
    public long STOP_AFTER = 0;

        /**
     * Final implementation of doWork() that checks and loads the input and optionally reference
     * sequence files and the runs the sublcass through the setup() acceptRead() and finish() steps.
     */
    @Override
    protected final int doWork() {
    	log.debug("Setting language-neutral locale");
    	java.util.Locale.setDefault(Locale.ROOT);
        try {
			makeItSo(INPUT, REFERENCE_SEQUENCE, ASSUME_SORTED, STOP_AFTER, Arrays.asList(this));
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e);
		}
        return 0;
    }

    public static void makeItSo(final File input,
                                final File referenceSequence,
                                final boolean assumeSorted,
                                final long stopAfter,
                                final Collection<ByReadNameSinglePassSamProgram> programs) throws FileNotFoundException {
        // Setup the standard inputs
        IOUtil.assertFileIsReadable(input);
        SamReader in = SamReaderFactory.makeDefault().referenceSequence(referenceSequence).open(input);
        // Optionally load up the reference sequence and double check sequence dictionaries
        final ReferenceLookup lookup;
        if (referenceSequence == null) {
        	lookup = null;
        } else {
            IOUtil.assertFileIsReadable(referenceSequence);
            lookup = new TwoBitBufferedReferenceSequenceFile(new IndexedFastaSequenceFile(referenceSequence));

            if (!in.getFileHeader().getSequenceDictionary().isEmpty()) {
                SequenceUtil.assertSequenceDictionariesEqual(in.getFileHeader().getSequenceDictionary(),
                		lookup.getSequenceDictionary());
            }
        }
        // Check on the sort order of the BAM file
        final SortOrder sort = in.getFileHeader().getSortOrder();
        if (sort != SortOrder.queryname) {
            if (assumeSorted) {
                log.warn("File reports sort order '" + sort + "', assuming it's queryname sorted anyway.");
            } else {
                throw new PicardException("File " + input.getAbsolutePath() + " should be queryname sorted but " +
                        "the header says the sort order is " + sort + ". If you believe the file " +
                        "to be queryname sorted you may pass ASSUME_SORTED=true");
            }
        }
        for (final ByReadNameSinglePassSamProgram program : programs) {
        	program.setReference(lookup);
            program.setup(in.getFileHeader(), input);
        }
        final ProgressLogger progress = new ProgressLogger(log);
        final SAMRecordIterator rawit = in.iterator();
        final CloseableIterator<SAMRecord> it = new AsyncBufferedIterator<SAMRecord>(rawit, "ByReadNameSinglePassSamProgram " + input.getName());
        try {
	        List<SAMRecord> currentRecords = new ArrayList<>();
	        String currentReadName = null;
	        while (it.hasNext()) {
	        	SAMRecord r = it.next();
	        	String readname = r.getReadName();
	        	// if read name we have to just treat it as a single read
	        	if (readname == null || !readname.equals(currentReadName)) {
	        		if (currentRecords.size() > 0) {
		        		for (final ByReadNameSinglePassSamProgram program : programs) {
		        			program.acceptFragment(currentRecords, lookup);
			            }
	        		}
	        		currentRecords.clear();
	        		currentReadName = readname;
	        		if (stopAfter > 0 && progress.getCount() >= stopAfter) {
		                break;
		            }
	        	}
	        	currentRecords.add(r);
	        	progress.record(r);
	        }
	        if (currentRecords.size() > 0) {
	        	for (final ByReadNameSinglePassSamProgram program : programs) {
	    			program.acceptFragment(currentRecords, lookup);
	            }
	        }
        } finally {
	        CloserUtil.close(it);
	        CloserUtil.close(rawit);
	        CloserUtil.close(in);
        }
        for (final ByReadNameSinglePassSamProgram program : programs) {
            program.finish();
        }
    }
    /** Should be implemented by subclasses to do one-time initialization work. */
    protected abstract void setup(final SAMFileHeader header, final File samFile);
    /**
     * Should be implemented by subclasses to accept SAMRecords one fragment at at time.
     * If a reference sequence file was supplied to the program it will be passed as 'ref'. Otherwise 'ref' may be null.
     */
    protected abstract void acceptFragment(final List<SAMRecord> records, ReferenceLookup lookup);
    /** Should be implemented by subclasses to do one-time finalization work. */
    protected abstract void finish();
    public void copyInput(ProcessStructuralVariantReadsCommandLineProgram to) {
    	CommandLineProgramHelper.copyInputs(this, to);
    	to.INPUT = INPUT;
    	to.OUTPUT = OUTPUT;
    	to.ASSUME_SORTED = ASSUME_SORTED;
    	to.STOP_AFTER = STOP_AFTER;
    }
}
