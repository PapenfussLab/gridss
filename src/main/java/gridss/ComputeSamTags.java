package gridss;

import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.sam.TemplateTagsIterator;
import au.edu.wehi.idsv.util.FileHelper;
import au.edu.wehi.idsv.util.ParallelTransformIterator;
import au.edu.wehi.idsv.util.UngroupingIterator;
import com.google.common.collect.Sets;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import gridss.cmdline.ReferenceCommandLineProgram;
import htsjdk.samtools.*;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import java.util.Set;
import java.util.concurrent.Executor;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

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
			//SAMTag.Q2.name(), // dropping Q2 to improve runtime performance and file size
			SAMTag.MC.name(),
			SAMTag.MQ.name());
	@Argument(doc = "Number of worker threads to spawn. Defaults to number of cores available."
			+ " Note that I/O threads are not included in this worker thread count so CPU usage can be higher than the number of worker thread.",
			shortName = "THREADS")
	public int WORKER_THREADS = Runtime.getRuntime().availableProcessors();
	@Override
	protected int doWork() {
		log.debug("Setting language-neutral locale");
    	java.util.Locale.setDefault(Locale.ROOT);
    	validateParameters();
    	SamReaderFactory readerFactory = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE);
    	SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
    	try {
			ensureDictionaryMatches(INPUT);
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
				ProgressLogger progress = new ProgressLogger(log);
				ExecutorService threadpool = Executors.newFixedThreadPool(WORKER_THREADS, new ThreadFactoryBuilder().setDaemon(false).setNameFormat("ComputeSamTags-%d").build());
    			try (SAMRecordIterator it = reader.iterator()) {
					Iterator<SAMRecord> asyncIt = transform(threadpool, Defaults.ASYNC_BUFFER_SIZE, it);
    				File tmpoutput = gridss.Defaults.OUTPUT_TO_TEMP_FILE ? FileSystemContext.getWorkingFileFor(OUTPUT, "gridss.tmp.ComputeSamTags.") : OUTPUT;
    				try (SAMFileWriter writer = writerFactory.makeSAMOrBAMWriter(header, true, tmpoutput)) {
						while (asyncIt.hasNext()) {
							SAMRecord r = asyncIt.next();
							writer.addAlignment(r);
							progress.record(r);
						}
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
	public Iterator<SAMRecord> transform(Executor threadpool, int batchSize, Iterator<SAMRecord> it) {
		final Set<String> tags = TAGS;
		final boolean softenHardClips = SOFTEN_HARD_CLIPS;
		final boolean fixMates = FIX_MATE_INFORMATION;
		final boolean fixDuplicates = FIX_DUPLICATE_FLAG;
		final boolean fixSA = FIX_SA;
		final boolean fixTruncated = FIX_MISSING_HARD_CLIP;
		final boolean recalculateSupplementary = RECALCULATE_SA_SUPPLEMENTARY;
		final ReferenceLookup reference = isReferenceRequired() ? getReference() : null;
		Iterator<List<SAMRecord>> groupByFragment = TemplateTagsIterator.withGrouping(it);
		ParallelTransformIterator<List<SAMRecord>, List<SAMRecord>> parallelIt = new ParallelTransformIterator<>(
				groupByFragment,
				r -> transform(
					r,
					tags,
					softenHardClips,
					fixMates,
					fixDuplicates,
					fixSA,
					fixTruncated,
					recalculateSupplementary,
					reference),
				batchSize,
				threadpool);
		return new UngroupingIterator(parallelIt);
	}
	private static List<SAMRecord> transform(
			List<SAMRecord> records,
			Set<String> tags,
			boolean softenHardClips,
			boolean fixMates,
			boolean fixDuplicates,
			boolean fixSA,
			boolean fixTruncated,
			boolean recalculateSupplementary,
			ReferenceLookup reference) {
		if ((tags.contains(SAMTag.NM.name()) || tags.contains(SAMTag.SA.name()))) {
			for (SAMRecord r : records) {
				SAMRecordUtil.ensureNmTag(reference, r);
			}
		}
		SAMRecordUtil.calculateTemplateTags(records, tags, softenHardClips, fixMates, fixDuplicates, fixSA, fixTruncated, recalculateSupplementary);
		return records;
	}
	protected boolean isReferenceRequired() {
		return TAGS.contains(SAMTag.NM.name()) ||
				TAGS.contains(SAMTag.SA.name()); // SA requires NM
	}
	protected void validateParameters() {
		if (isReferenceRequired()) {
			IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
		}
    	IOUtil.assertFileIsReadable(INPUT);
    	IOUtil.assertFileIsWritable(OUTPUT);
	}
	@Override
	protected String[] customCommandLineValidation() {
		if (isReferenceRequired()) {
			if (REFERENCE_SEQUENCE == null) {
				return new String[] { "Must have a non-null REFERENCE_SEQUENCE"};
			}
			if (!REFERENCE_SEQUENCE.exists()) {
				return new String[] { "Missing REFERENCE_SEQUENCE " + REFERENCE_SEQUENCE};
			}
        }
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
						return new String[] { msg } ;
				}
			} catch (IllegalArgumentException e) {
				// ignore
			}
		}
		return super.customCommandLineValidation();
	}
	public static void main(String[] argv) {
        System.exit(new ComputeSamTags().instanceMain(argv));
    }
}
