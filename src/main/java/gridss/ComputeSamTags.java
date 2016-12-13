package gridss;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.Locale;
import java.util.Set;

import com.google.common.collect.Sets;

import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.sam.NmTagIterator;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.sam.TemplateTagsIterator;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import au.edu.wehi.idsv.util.FileHelper;
import gridss.cmdline.ReferenceCommandLineProgram;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(
        usage = "Populates computed SAM tags. "
        		+ "The NM tags requires the reference genome to be specified. "
        		+ "Tags requiring information from mate reads, or alternative, secondary, or chimeric alignments require"
        		+ " all reads from the same fragment/template to be sorted together."
        		+ " This can be achieved by queryname sorting of the input file, although the raw output from most aligners"
        		+ " also fulfills this criteria.",
        usageShort = "Populates computed SAM tags."
)
public class ComputeSamTags extends ReferenceCommandLineProgram {
	private static final Log log = Log.getInstance(ComputeSamTags.class);
	@Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input BAM file grouped by read name.")
    public File INPUT;
	@Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Annotated BAM file.")
    public File OUTPUT;
	@Option(shortName=StandardOptionDefinitions.ASSUME_SORTED_SHORT_NAME, doc="Assume that all records with the same read name are consecutive. "
			+ "Incorrect tags will be written if this is not the case.", optional=true)
    public boolean ASSUME_SORTED = false;
	@Option(doc="Convert hard clips to soft clips if the entire read sequence for the read is available in another record.", optional=true)
	public boolean SOFTEN_HARD_CLIPS = true;
	@Option(doc="Fixes missing mate information. Unlike Picard tools FixMateInformation, reads for which no mate can be found"
			+ " are converted to unpaired reads.", optional=true)
	public boolean FIX_MATE_INFORMATION = true;
	@Option(shortName="T", doc="Tags to calculate")
	public Set<SAMTag> TAGS = Sets.newHashSet(
			SAMTag.NM,
			SAMTag.SA,
			SAMTag.Q2,
			SAMTag.R2,
			SAMTag.MC,
			SAMTag.MQ);
	@Override
	protected int doWork() {
		log.debug("Setting language-neutral locale");
    	java.util.Locale.setDefault(Locale.ROOT);
    	validateParameters();
    	SamReaderFactory readerFactory = SamReaderFactory.make();
    	SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
    	try {
    		try (SamReader reader = readerFactory.open(INPUT)) {
    			SAMFileHeader header = reader.getFileHeader();
    			if (!ASSUME_SORTED) {
    				if (header.getSortOrder() != SortOrder.queryname) {
    					log.error("INPUT is not sorted by queryname. "
    							+ "ComputeSamTags requires that reads with the same name be sorted together. "
    							+ "If the input file satisfies this constraint (the output from many aligners do),"
    							+ " this check can be disabled with the ASSUME_SORTED option.");
    					return -1;
    				}
    			}
    			try (SAMRecordIterator it = reader.iterator()) {
    				File tmpoutput = FileSystemContext.getWorkingFileFor(OUTPUT, "gridss.tmp.ComputeSamTags.");
    				try (SAMFileWriter writer = writerFactory.makeSAMOrBAMWriter(header, true, tmpoutput)) {
    					compute(it, writer, getReference(), TAGS, SOFTEN_HARD_CLIPS, FIX_MATE_INFORMATION, INPUT.getName() + "-");
    				}
    				FileHelper.move(tmpoutput, OUTPUT, true);
    			}
    		}
		} catch (IOException e) {
			log.error(e);
			return -1;
		}
    	return 0;
	}
	public static void compute(Iterator<SAMRecord> rawit, SAMFileWriter writer, ReferenceLookup reference, Set<SAMTag> tags, boolean softenHardClips, boolean fixMates, String threadprefix) throws IOException {
		ProgressLogger progress = new ProgressLogger(log);
		try (CloseableIterator<SAMRecord> aysncit = new AsyncBufferedIterator<SAMRecord>(rawit, threadprefix + "raw")) {
			Iterator<SAMRecord> it = aysncit;
			if (tags.contains(SAMTag.NM) || tags.contains(SAMTag.SA)) {
				it = new AsyncBufferedIterator<SAMRecord>(it, threadprefix + "nm");
				it = new NmTagIterator(it, reference);
			}
			if (!Sets.intersection(tags, SAMRecordUtil.TEMPLATE_TAGS).isEmpty() || softenHardClips) {
				it = new TemplateTagsIterator(it, softenHardClips, fixMates, tags);
				it = new AsyncBufferedIterator<SAMRecord>(it, threadprefix + "tags");
			}
			while (it.hasNext()) {
				SAMRecord r = it.next();
				writer.addAlignment(r);
				progress.record(r);
			}
		}
	}
	private boolean isReferenceRequired() {
		return TAGS.contains(SAMTag.NM) ||
				TAGS.contains(SAMTag.SA); // SA requires NM
	}
	private void validateParameters() {
		if (isReferenceRequired()) {
			IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
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
    		case MC:
    		case MQ:
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
	@Override
	protected String[] customCommandLineValidation() {
		if (isReferenceRequired() && REFERENCE_SEQUENCE == null) {
            return new String[]{"Must have a non-null REFERENCE_SEQUENCE"};
        }
		return super.customCommandLineValidation();
	}
	public static void main(String[] argv) {
        System.exit(new ComputeSamTags().instanceMain(argv));
    }
}
