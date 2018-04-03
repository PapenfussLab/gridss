package gridss;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;

import com.google.common.collect.Iterators;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.DirectedBreakpoint;
import au.edu.wehi.idsv.IndelEvidence;
import au.edu.wehi.idsv.NonReferenceReadPair;
import au.edu.wehi.idsv.ProgressLoggingSAMRecordIterator;
import au.edu.wehi.idsv.SingleReadEvidence;
import au.edu.wehi.idsv.SoftClipEvidence;
import au.edu.wehi.idsv.SplitReadEvidence;
import au.edu.wehi.idsv.StringEvidenceIdentifierGenerator;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import au.edu.wehi.idsv.util.MathUtil;
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
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(
		summary = "Converts split reads and indel-containing reads to BEDPE notation.",
		oneLineSummary = "Converts split reads and indel-containing reads to BEDPE notation.",
        programGroup = gridss.cmdline.programgroups.DataConversion.class)
public class ReadsToBedpe extends CommandLineProgram {
	private static final Log log = Log.getInstance(ReadsToBedpe.class);
    @Argument(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input file", optional=false)
    public File INPUT;
    @Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output BEDPE", optional=false)
    public File OUTPUT;
    @Argument(doc="Minimum event size", optional=true)
    public int MIN_SIZE = 50;
    @Argument(doc="Minimum mapq of read alignments", optional=true)
    public int MIN_MAPQ = 0;
    @Argument(doc="Write split reads", optional=true)
    public boolean SPLIT_READS = true;
    @Argument(doc="Write indel reads", optional=true)
    public boolean INDELS = true;
    @Argument(doc="Value to write to the BEDPE name field. Note that the unique identifier includes the read CIGAR so can be very long for long read sequencing technologies.", optional=true)
    public Name NAME = Name.ReadName;
    /**
     * Value to populate the name field
     * @author Daniel Cameron
     *
     */
    public enum Name {
    	None,
    	ReadName,
    	UniqueIdentifier,
    }
    @Override
	protected int doWork() {
		log.debug("Setting language-neutral locale");
    	java.util.Locale.setDefault(Locale.ROOT);
    	validateParameters();
    	SamReaderFactory readerFactory = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE);
    	try {
    		try (SamReader reader = readerFactory.open(INPUT)) {
    			SAMFileHeader header = reader.getFileHeader();
    			SAMSequenceDictionary dict = header.getSequenceDictionary();
    			//log.info(String.format("Using %d worker threads", WORKER_THREADS));
        		//ExecutorService threadpool = Executors.newFixedThreadPool(WORKER_THREADS, new ThreadFactoryBuilder().setDaemon(false).setNameFormat("Worker-%d").build());
    			try (CloseableIterator<SAMRecord> rawit = new AsyncBufferedIterator<SAMRecord>(reader.iterator(), 3, 64)) {
    				ProgressLoggingSAMRecordIterator logit = new ProgressLoggingSAMRecordIterator(rawit, new ProgressLogger(log));
    				//ParallelTransformIterator<SAMRecord, List<String>> it = new ParallelTransformIterator<>(logit, r -> asBedPe(dict, r), 16 + 2 * WORKER_THREADS, threadpool);
    				Iterator<List<String>> it = Iterators.transform(logit, r -> asBedPe(dict, r));
    				
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
		for (IndelEvidence ie : IndelEvidence.create(null, MIN_SIZE, record)) {
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
		sb.append(Integer.toString(bp.start - 1));
		sb.append('	');
		sb.append(Integer.toString(bp.end));
		sb.append('	');
		sb.append(dict.getSequence(bp.referenceIndex2).getSequenceName());
		sb.append('	');
		sb.append(Integer.toString(bp.start2 - 1));
		sb.append('	');
		sb.append(Integer.toString(bp.end2));
		sb.append('	');
		switch (NAME) {
		case ReadName:
			if (e instanceof SingleReadEvidence) {
				sb.append(((SingleReadEvidence)e).getSAMRecord().getReadName());
			} else {
				sb.append('.');
			}
			break;
		case UniqueIdentifier:
			if (e instanceof NonReferenceReadPair) {
				sb.append(new StringEvidenceIdentifierGenerator().getEvidenceID((NonReferenceReadPair)e));
			} else if (e instanceof SoftClipEvidence) {
				sb.append(new StringEvidenceIdentifierGenerator().getEvidenceID((SoftClipEvidence)e));
			} else if (e instanceof SplitReadEvidence) {
				sb.append(new StringEvidenceIdentifierGenerator().getEvidenceID((SplitReadEvidence)e));
			} else if (e instanceof IndelEvidence) {
				sb.append(new StringEvidenceIdentifierGenerator().getEvidenceID((IndelEvidence)e));
			} else {
				sb.append('.');
			}
			break;
		default:
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
