package gridss;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.IOException;
import java.util.Locale;

import org.apache.commons.lang3.NotImplementedException;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import au.edu.wehi.idsv.DirectedBreakpoint;
import au.edu.wehi.idsv.IndelEvidence;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;

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
    			//SAMFileHeader header = reader.getFileHeader();
    			try (CloseableIterator<SAMRecord> it = new AsyncBufferedIterator<SAMRecord>(reader.iterator(), 2, 300)) {
    				process(it.next());
    			}
    		}
		} catch (IOException e) {
			log.error(e);
			return -1;
		}
    	return 0;
	}
	private void process(SAMRecord record) {
		// Split read
		// indel
		for (IndelEvidence ie : IndelEvidence.create(null, record)) {
			writeBedPe(ie, "indel");
		}
	}
	private void writeBedPe(DirectedBreakpoint bp, String source) {
		throw new NotImplementedException();
	}
	private void validateParameters() {
    	IOUtil.assertFileIsReadable(INPUT);
    	IOUtil.assertFileIsWritable(OUTPUT);
	}
	public static void main(String[] argv) {
        System.exit(new ReadsToBedpe().instanceMain(argv));
    }
}
