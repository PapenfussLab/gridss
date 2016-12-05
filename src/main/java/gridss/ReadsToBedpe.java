package gridss;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.Locale;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.DirectedBreakpoint;
import au.edu.wehi.idsv.IndelEvidence;
import au.edu.wehi.idsv.SplitReadEvidence;
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
    //@Option(shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="Reference genome", optional=true)
    //public File REFERENCE = null;
    @Option(doc="Minimum event size", optional=true)
    public int MIN_SIZE = 50;
    @Option(doc="Minimum mapq of read alignments", optional=true)
    public int MIN_MAPQ = 0;
    @Option(doc="Write split reads", optional=true)
    public boolean SPLIT_READS = true;
    @Option(doc="Write indel reads", optional=true)
    public boolean INDELS = true;
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
    			try (CloseableIterator<SAMRecord> it = new AsyncBufferedIterator<SAMRecord>(reader.iterator(), 2, 300)) {
    				int i= 0;
    				try (BufferedWriter writer = new BufferedWriter(new FileWriter(OUTPUT))) {
    					while (it.hasNext()) {
    						process(writer, dict, it.next());
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
	private void process(Writer writer, SAMSequenceDictionary dict, SAMRecord record) throws IOException {
		// Split read
		for (SplitReadEvidence sre : SplitReadEvidence.create(null, record)) {
			if (sre.getBreakendSummary().isLowBreakend()) {
				writeBedPe(writer, dict, sre);
			}
		}
		// indel
		for (IndelEvidence ie : IndelEvidence.create(null, record)) {
			if (ie.getBreakendSummary().direction == BreakendDirection.Forward) {
				writeBedPe(writer, dict, ie);
			}
		}
	}
	private void writeBedPe(Writer writer, SAMSequenceDictionary dict, SplitReadEvidence e) throws IOException {
		writeBedPe(writer, dict, e, (int)MathUtil.phredOr(e.getLocalMapq(), e.getRemoteMapq()), "splitread");
	}
	private void writeBedPe(Writer writer, SAMSequenceDictionary dict, IndelEvidence e) throws IOException {
		writeBedPe(writer, dict, e, e.getLocalMapq(), "indel");
	}
	private void writeBedPe(Writer writer, SAMSequenceDictionary dict, DirectedBreakpoint e, int mapq, String source) throws IOException {
		Integer size = e.getBreakendSummary().getEventSize();
		if (size == null || Math.abs(size) < MIN_SIZE) return;
		if (e.getLocalMapq() < MIN_MAPQ || e.getRemoteMapq() < MIN_MAPQ) return;
		
		BreakpointSummary bp = e.getBreakendSummary();
		writer.write(dict.getSequence(bp.referenceIndex).getSequenceName());
		writer.write('	');
		writer.write(Integer.toString(bp.start));
		writer.write('	');
		writer.write(Integer.toString(bp.end));
		writer.write('	');
		writer.write(dict.getSequence(bp.referenceIndex2).getSequenceName());
		writer.write('	');
		writer.write(Integer.toString(bp.start2));
		writer.write('	');
		writer.write(Integer.toString(bp.end2));
		writer.write('	');
		writer.write(e.getEvidenceID());
		writer.write('	');
		writer.write(Integer.toString(mapq));
		writer.write('	');
		writer.write(bp.direction == BreakendDirection.Forward ? '+' : '-');
		writer.write('	');
		writer.write(bp.direction2 == BreakendDirection.Forward ? '+' : '-');
		writer.write('	');
		writer.write(source);
		writer.write('\n');
	}
	private void validateParameters() {
    	IOUtil.assertFileIsReadable(INPUT);
    	IOUtil.assertFileIsWritable(OUTPUT);
	}
	public static void main(String[] argv) {
        System.exit(new ReadsToBedpe().instanceMain(argv));
    }
}
