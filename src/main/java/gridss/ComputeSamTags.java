package gridss;

import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.sam.NmTagIterator;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.sam.TemplateTagsIterator;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import au.edu.wehi.idsv.util.FileHelper;
import com.google.common.collect.Sets;
import gridss.cmdline.ReferenceCommandLineProgram;
import htsjdk.samtools.*;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.Locale;
import java.util.Set;

@CommandLineProgramProperties(
		summary = "Populates computed SAM tags. "
        		+ "The NM tags requires the reference genome to be specified. "
        		+ "Tags requiring information from mate reads, or alternative, secondary, or chimeric alignments require"
        		+ " all reads from the same fragment/template to be sorted together."
        		+ " This can be achieved by queryname sorting of the input file, although the raw output from most aligners"
        		+ " also fulfills this criteria.",
        oneLineSummary = "Populates computed SAM tags.",
        programGroup = picard.cmdline.programgroups.ReadDataManipulationProgramGroup.class
)
public class ComputeSamTags extends ReferenceCommandLineProgram {
	private static final Log log = Log.getInstance(ComputeSamTags.class);
	@Argument(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input BAM file grouped by read name.")
    public File INPUT;
	@Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Annotated BAM file.")
    public File OUTPUT;
	@Argument(shortName=StandardOptionDefinitions.ASSUME_SORTED_SHORT_NAME, doc="Assume that all records with the same read name are consecutive. "
			+ "Incorrect tags will be written if this is not the case.", optional=true)
    public boolean ASSUME_SORTED = false;
	@Argument(doc="Convert hard clips to soft clips if the entire read sequence for the read is available in another record. "
			+ "If no base information can be found, N bases with 0 base quality are substituted.", optional=true)
	public boolean SOFTEN_HARD_CLIPS = true;
	@Argument(doc="Fixes missing mate information. Unlike Picard tools FixMateInformation, reads for which no mate can be found"
			+ " are converted to unpaired reads.", optional=true)
	public boolean FIX_MATE_INFORMATION = true;
	@Argument(doc="Sets the duplicate flag if any alignment in the read pair is flagged as a duplicate. Many duplicate marking tools do not correctly mark all supplementary alignments.", optional=true)
	public boolean FIX_DUPLICATE_FLAG = true;
	@Argument(doc="Fixes the SA tag to match the read alignments. Useful for programs such as GATK indel realignment do not update the SA tag when adjusting read alignments.", optional=true)
	public boolean FIX_SA = true;
	@Argument(doc="Adds hard clipping CIGAR elements to truncated alignments. Useful for programs such as GATK indel realignment that strip hard clips. Assumes all alignments form part of the split read thus does not support secondary alignments.", optional=true)
	public boolean FIX_MISSING_HARD_CLIP = true;
	@Argument(doc="Recalculates the supplementary flag based on the SA tag. The supplementary flag should be set on all split read alignments except one.", optional=true)
	public boolean RECALCULATE_SA_SUPPLEMENTARY = true;
	@Argument(shortName="T", doc="Tags to calculate")
	public Set<String> TAGS = Sets.newHashSet(
			SAMTag.NM.name(),
			SAMTag.SA.name(),
			SAMTag.R2.name(),
			//SAMTag.Q2.name(), // dropping Q2 to improve performance
			SAMTag.MC.name(),
			SAMTag.MQ.name());
	@Override
	protected int doWork() {
		log.debug("Setting language-neutral locale");
    	java.util.Locale.setDefault(Locale.ROOT);
    	validateParameters();
    	SamReaderFactory readerFactory = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE);
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
    			} else {
    				// strip sort order header so we can use samtools sorted files
    				header.setSortOrder(SortOrder.unsorted);
    			}
    			try (SAMRecordIterator it = reader.iterator()) {
    				File tmpoutput = gridss.Defaults.OUTPUT_TO_TEMP_FILE ? FileSystemContext.getWorkingFileFor(OUTPUT, "gridss.tmp.ComputeSamTags.") : OUTPUT;
    				try (SAMFileWriter writer = writerFactory.makeSAMOrBAMWriter(header, true, tmpoutput)) {
    					compute(it, writer, getReference(), TAGS, SOFTEN_HARD_CLIPS, FIX_MATE_INFORMATION, FIX_DUPLICATE_FLAG, FIX_SA, FIX_MISSING_HARD_CLIP, RECALCULATE_SA_SUPPLEMENTARY, INPUT.getName() + "-");
    				}
    				if (tmpoutput != OUTPUT) {
    					FileHelper.move(tmpoutput, OUTPUT, true);
    				}
    			}
    		}
		} catch (IOException e) {
			log.error(e);
			return -1;
		}
    	return 0;
	}
	public static void compute(Iterator<SAMRecord> rawit, SAMFileWriter writer, ReferenceLookup reference, Set<String> tags,
			boolean softenHardClips,
			boolean fixMates,
			boolean fixDuplicates,
		    boolean fixSA,
		    boolean fixTruncated,
			boolean recalculateSupplementary,
			String threadprefix) throws IOException {
		ProgressLogger progress = new ProgressLogger(log);
		try (CloseableIterator<SAMRecord> aysncit = new AsyncBufferedIterator<SAMRecord>(rawit, threadprefix + "raw")) {
			Iterator<SAMRecord> it = aysncit;
			if (tags.contains(SAMTag.NM.name()) || tags.contains(SAMTag.SA.name())) {
				it = new AsyncBufferedIterator<SAMRecord>(it, threadprefix + "nm");
				it = new NmTagIterator(it, reference);
			}
			if (!Sets.intersection(tags, SAMRecordUtil.TEMPLATE_TAGS).isEmpty() || softenHardClips) {
				it = new TemplateTagsIterator(it, softenHardClips, fixMates, fixDuplicates, fixSA, fixTruncated, recalculateSupplementary, tags);
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
		return TAGS.contains(SAMTag.NM.name()) ||
				TAGS.contains(SAMTag.SA.name()); // SA requires NM
	}
	private void validateParameters() {
		if (isReferenceRequired()) {
			IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
		}
    	IOUtil.assertFileIsReadable(INPUT);
    	IOUtil.assertFileIsWritable(OUTPUT);
    	for (String t : TAGS) {
    		try {
    			SAMTag tag = SAMTag.valueOf(t);
    			switch (tag) {
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
    		} catch (IllegalArgumentException e) {
    			// ignore
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
