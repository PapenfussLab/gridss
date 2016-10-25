package gridss;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Locale;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.GenomicProcessingContext;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;

@CommandLineProgramProperties(
		usage = "Converts split reads and indel-containing reads to BEDPE notation.", usageShort = "")
public class ReadsToBedpe extends CommandLineProgram {
	private static final Log log = Log.getInstance(ComputeSamTags.class);
    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input file", optional=false)
    public File INPUT;
    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Input BEDPE", optional=false)
    public File OUTPUT;
    @Option(shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="Reference genome", optional=true)
    public File REFERENCE = null;
    @Option(doc="Minimum event size", optional=true)
    public int MIN_SIZE = 50;
    @Option(doc="Write split reads", optional=true)
    public boolean SPLIT_READS = true;
    @Option(doc="Write indel reads", optional=true)
    public boolean INDELS = true;
    private ReferenceSequenceFile reference;
    private ReferenceSequenceFile getReference() throws FileNotFoundException {
    	if (reference == null) {
    		reference = new IndexedFastaSequenceFile(REFERENCE);
    	}
    	return reference;
    }
    @Override
	protected int doWork() {
		log.debug("Setting language-neutral locale");
    	java.util.Locale.setDefault(Locale.ROOT);
    	validateParameters();
    	GenomicProcessingContext pc = new GenomicProcessingContext(new FileSystemContext(TMP_DIR.get(0), MAX_RECORDS_IN_RAM), REFERENCE, false, null);
    	SamReaderFactory readerFactory = SamReaderFactory.make();
    	try {
    		try (SamReader reader = readerFactory.open(INPUT)) {
    			SAMFileHeader header = reader.getFileHeader();
    			try (CloseableIterator<SAMRecord> it = new AsyncBufferedIterator<SAMRecord>(reader.iterator(), 2, 128)) {
    				
    			}
    		}
		} catch (IOException e) {
			log.error(e);
			return -1;
		}
    	return 0;
	}
	private void validateParameters() {
    	IOUtil.assertFileIsReadable(INPUT);
    	IOUtil.assertFileIsWritable(OUTPUT);
	}
	public static void main(String[] argv) {
        System.exit(new ReadsToBedpe().instanceMain(argv));
    }
}
