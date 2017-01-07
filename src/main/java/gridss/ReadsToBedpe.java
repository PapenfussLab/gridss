package gridss;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import com.google.common.util.concurrent.ThreadFactoryBuilder;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.DirectedBreakpoint;
import au.edu.wehi.idsv.IndelEvidence;
import au.edu.wehi.idsv.ProgressLoggingSAMRecordIterator;
import au.edu.wehi.idsv.SplitReadEvidence;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import au.edu.wehi.idsv.util.MathUtil;
import au.edu.wehi.idsv.util.ParallelTransformIterator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(
		usage = "Converts split reads and indel-containing reads to BEDPE notation.", usageShort = "")
public class ReadsToBedpe extends CommandLineProgram {
	private static final Log log = Log.getInstance(ComputeSamTags.class);
    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input file", optional=false)
    public File INPUT;
    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output BEDPE", optional=false)
    public File OUTPUT;
    @Option(doc="Minimum event size", optional=true)
    public int MIN_SIZE = 50;
    @Option(doc="Minimum mapq of read alignments", optional=true)
    public int MIN_MAPQ = 0;
    @Option(doc="Write split reads", optional=true)
    public boolean SPLIT_READS = true;
    @Option(doc="Write indel reads", optional=true)
    public boolean INDELS = true;
    @Option(doc="Include a unique identifier for each breakpoint supported by each read. "
    		+ "Note that this identified can be quite long for long read sequencing technologies.", optional=true)
    public boolean UNIQUE_IDENTIFIER= true;
    @Option(doc="Number of worker threads to spawn. Defaults to number of cores available."
			+ " Note that I/O threads are not included in this worker thread count so CPU usage can be higher than the number of worker thread.",
    		shortName="THREADS")
    public int WORKER_THREADS = Runtime.getRuntime().availableProcessors();
    @Override
    protected String[] customCommandLineValidation() {
    	String[] val = super.customCommandLineValidation();
    	if (val == null || val.length == 0) {
	    	if (WORKER_THREADS <= 0) {
	    		return new String[] { "WORKER_THREADS must be at least 1." }; 
	    	}
    	}
    	return val;
    }
    @Override
	protected int doWork() {
		log.debug("Setting language-neutral locale");
    	java.util.Locale.setDefault(Locale.ROOT);
    	validateParameters();
    	SamReaderFactory readerFactory = SamReaderFactory.make();
    	try {
    		try (SamReader reader = readerFactory.open(INPUT)) {
    			SAMFileHeader header = reader.getFileHeader();
    			SAMSequenceDictionary dict = header.getSequenceDictionary();
    			log.info(String.format("Using %d worker threads", WORKER_THREADS));
        		ExecutorService threadpool = Executors.newFixedThreadPool(WORKER_THREADS, new ThreadFactoryBuilder().setDaemon(false).setNameFormat("Worker-%d").build());
    			try (CloseableIterator<SAMRecord> rawit = new AsyncBufferedIterator<SAMRecord>(reader.iterator(), 3, 64)) {
    				ProgressLoggingSAMRecordIterator logit = new ProgressLoggingSAMRecordIterator(rawit, new ProgressLogger(log));
    				ParallelTransformIterator<SAMRecord, List<String>> it = new ParallelTransformIterator<>(logit, r -> asBedPe(dict, r), 16 + 2 * WORKER_THREADS, threadpool);
    				int i = 0;
    				try (BufferedWriter writer = new BufferedWriter(new FileWriter(OUTPUT))) {
    					while (it.hasNext()) {
    						for (String line : it.next()) {
    							if (line != null) {
		    						writer.write(line);
		    						writer.write('\n');
    							}
    						}
    						i++;
    					}
    					if (i % 1000 == 0) {
    						writer.flush();
    					}
    				}
    			}
    		}
		} catch (IOException e) {
			log.error(e);
			return -1;
		}
    	return 0;
	}
	private List<String> asBedPe(SAMSequenceDictionary dict, SAMRecord record) {
		List<String> result = new ArrayList<>();
		// Split read
		for (SplitReadEvidence sre : SplitReadEvidence.create(null, record)) {
			if (sre.getBreakendSummary().isLowBreakend()) {
				result.add(asBedPe(dict, sre));
			}
		}
		// indel
		for (IndelEvidence ie : IndelEvidence.create(null, record)) {
			if (ie.getBreakendSummary().direction == BreakendDirection.Forward) {
				result.add(asBedPe(dict, ie));
			}
		}
		return result;
	}
	private String asBedPe(SAMSequenceDictionary dict, SplitReadEvidence e) {
		return asBedPe(dict, e, (int)MathUtil.phredOr(e.getLocalMapq(), e.getRemoteMapq()), "splitread");
	}
	private String asBedPe(SAMSequenceDictionary dict, IndelEvidence e) {
		return asBedPe(dict, e, e.getLocalMapq(), "indel");
	}
	private String asBedPe(SAMSequenceDictionary dict, DirectedBreakpoint e, int mapq, String source) {
		StringBuilder sb = new StringBuilder();
		Integer size = e.getBreakendSummary().getEventSize();
		if (size == null || Math.abs(size) < MIN_SIZE) return null;
		if (e.getLocalMapq() < MIN_MAPQ || e.getRemoteMapq() < MIN_MAPQ) return null;
		
		BreakpointSummary bp = e.getBreakendSummary();
		sb.append(dict.getSequence(bp.referenceIndex).getSequenceName());
		sb.append('	');
		sb.append(Integer.toString(bp.start));
		sb.append('	');
		sb.append(Integer.toString(bp.end));
		sb.append('	');
		sb.append(dict.getSequence(bp.referenceIndex2).getSequenceName());
		sb.append('	');
		sb.append(Integer.toString(bp.start2));
		sb.append('	');
		sb.append(Integer.toString(bp.end2));
		sb.append('	');
		if (UNIQUE_IDENTIFIER) {
			sb.append(e.getEvidenceID());
		} else {
			sb.append('.');
		}
		sb.append('	');
		sb.append(Integer.toString(mapq));
		sb.append('	');
		sb.append(bp.direction == BreakendDirection.Forward ? '+' : '-');
		sb.append('	');
		sb.append(bp.direction2 == BreakendDirection.Forward ? '+' : '-');
		sb.append('	');
		sb.append(source);
		sb.append('	');
		sb.append(e.getUntemplatedSequence().length());
		return sb.toString();
	}
	private void validateParameters() {
    	IOUtil.assertFileIsReadable(INPUT);
    	IOUtil.assertFileIsWritable(OUTPUT);
	}
	public static void main(String[] argv) {
        System.exit(new ReadsToBedpe().instanceMain(argv));
    }
}
