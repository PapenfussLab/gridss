package gridss;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.Locale;
import java.util.Set;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import au.edu.wehi.idsv.sam.NmTagIterator;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.sam.TemplateTagsIterator;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;

import com.google.common.collect.Sets;

@CommandLineProgramProperties(
        usage = "Populates computed SAM tags. "
        		+ "The NM tags requires the reference genome to be specified. "
        		+ "Tags requiring information from mate reads, or alternative, secondary, or chimeric alignments require"
        		+ " all reads from the same fragment/template to be sorted together."
        		+ " This can be achieved by queryname sorting of the input file, although the raw output from most aligners"
        		+ " also fulfills this criteria.",
        usageShort = "Populates computed SAM tags."
)
public class ComputeSamTags extends CommandLineProgram {
	private static final Log log = Log.getInstance(ComputeSamTags.class);
	@Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input BAM file grouped by read name.")
    public File INPUT;
	@Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Annotated BAM file.")
    public File OUTPUT;
	@Option(shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="Reference genome", optional=true)
    public File REFERENCE = null;
	@Option(shortName="T", doc="Tags to populate")
	public Set<SAMTag> TAGS = Sets.newHashSet(
			SAMTag.NM,
			SAMTag.Q2,
			SAMTag.R2,
			SAMTag.SA);
	@Override
	protected int doWork() {
		log.debug("Setting language-neutral locale");
    	java.util.Locale.setDefault(Locale.ROOT);
    	validateParameters();
    	SamReaderFactory readerFactory = SamReaderFactory.make().enable(htsjdk.samtools.SamReaderFactory.Option.EAGERLY_DECODE);
    	SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
    	try {
    		ReferenceSequenceFile reference = null;
    		if (isReferenceRequired()) {
    			reference = new IndexedFastaSequenceFile(REFERENCE);
    		}
			compute(INPUT, OUTPUT, reference, readerFactory, writerFactory, TAGS);
		} catch (IOException e) {
			log.error(e);
			return -1;
		}
    	return 0;
	}
	public static void compute(File input, File output, ReferenceSequenceFile reference, SamReaderFactory readerFactory, SAMFileWriterFactory writerFactory, Set<SAMTag> tags) throws IOException {
		try (SamReader reader = readerFactory.open(input)) {
			try (SAMRecordIterator it = reader.iterator()) {
				try (SAMFileWriter writer = writerFactory.makeSAMOrBAMWriter(reader.getFileHeader(), true, output)) {
					compute(it, writer, reference, tags);
				}
			}
		}
	}
	public static void compute(Iterator<SAMRecord> rawit, SAMFileWriter writer, ReferenceSequenceFile reference, Set<SAMTag> tags) throws IOException {
		try (CloseableIterator<SAMRecord> aysncit = new AsyncBufferedIterator<SAMRecord>(rawit, "raw records", 2, 300)) {
			Iterator<SAMRecord> it = aysncit;
			if (tags.contains(SAMTag.NM) || tags.contains(SAMTag.SA)) {
				it = new NmTagIterator(it, reference);
			}
			if (!Sets.intersection(tags, SAMRecordUtil.TEMPLATE_TAGS).isEmpty()) {
				it = new TemplateTagsIterator(it);
			}
			it = new AsyncBufferedIterator<SAMRecord>(it, "annotation", 2, 300);
			while (it.hasNext()) {
				writer.addAlignment(it.next());
			}
		}
	}
	private boolean isReferenceRequired() {
		return TAGS.contains(SAMTag.NM) ||
				TAGS.contains(SAMTag.SA); // SA requires NM
	}
	private void validateParameters() {
		if (isReferenceRequired()) {
			IOUtil.assertFileIsReadable(REFERENCE);
		}
    	IOUtil.assertFileIsReadable(INPUT);
    	IOUtil.assertFileIsWritable(OUTPUT);
    	for (SAMTag t : TAGS) {
    		switch (t) {
    		case CC:
    		case CP:
    		case FI:
    		case HI:
    		case IH:
    		case NM:
    		case Q2:
    		case R2:
    		case SA:
    		case TC:
    			break;
			default:
				String msg = String.format("%s is not a predefined standard SAM tag able to be computed with no additional information.", t); 
				log.error(msg);
				throw new RuntimeException(msg);
    		}
    	}
	}
	public static void main(String[] argv) {
        System.exit(new ComputeSamTags().instanceMain(argv));
    }
}
