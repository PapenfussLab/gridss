package gridss;

import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.SAMRecordChangeTracker;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import au.edu.wehi.idsv.util.FileHelper;
import au.edu.wehi.idsv.util.ParallelTransformIterator;
import au.edu.wehi.idsv.util.UngroupingIterator;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import gridss.cmdline.ReferenceCommandLineProgram;
import htsjdk.samtools.*;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.Executor;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import static au.edu.wehi.idsv.sam.SAMRecordUtil.ensureConsistentDuplicateFlag;

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
	@Argument(doc="Adds hard clipping CIGAR elements to truncated alignments." +
			"Useful for programs such as GATK indel realignment that strip hard clips." +
			"Hard clipping position is inferred based on the sequence of unclipped/soft clipped records for that reads.", optional=true)
	public boolean FIX_MISSING_HARD_CLIP = true;
	@Argument(doc="Recalculates the supplementary flag based on the SA tag. " +
			"The supplementary flag should be set on all split read alignments except one. " +
			"WARNING: this option treats secondary alignments as supplementary due to many pipelines still using 'bwa mem -M'.", optional=true)
	public boolean RECALCULATE_SUPPLEMENTARY = true;
	@Argument(doc="Alignment overlap threshold (in bp) for determining whether a read is a valid chimeric record or should be filtered out." +
			"Corrects issues such as" +
			"minimap2 reporting supplementary alignments fully overlapping the primary" +
			"and bwa reporting ALT contig alignments violating the SAM specifications (https://github.com/lh3/bwa/issues/282).", optional=true)
	public int SUPPLEMENTARY_ALIGNMENT_OVERLAP_THRESHOLD = 25;
	@Argument(doc="Fix CIGARs that have I or D closer to the end of the read than an aligned base." +
			"Since bwa does not follow the SAM specifications regarding alignment start position on these records, the position with the lowest alignment edit distance is chosen", optional=true)
	public boolean FIX_TERMINAL_CIGAR_INDEL = true;
	@Argument(doc="Outputs a tsv containing an overview of the changes made.", optional=true)
	public File MODIFICATION_SUMMARY_FILE = null;
	@Argument(shortName="T", doc="Tags to recalculate. SA recalculation is useful for programs such as GATK indel realignment do not update the SA tag when adjusting read alignments.")
	public Set<String> TAGS = Sets.newHashSet(
			SAMTag.NM.name(),
			SAMTag.SA.name(),
			SAMTag.R2.name(),
			//SAMTag.Q2.name(), // dropping Q2 to improve runtime performance and file size
			SAMTag.MC.name(),
			SAMTag.MQ.name());
	@Argument(doc="Tags to remove.")
	public List<String> REMOVE_TAGS = null;
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
    	if (REMOVE_TAGS == null) {
    		REMOVE_TAGS = Collections.emptyList();
		}
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
				SAMRecordChangeTracker tracker = null;
    			if (MODIFICATION_SUMMARY_FILE != null) {
					tracker = new SAMRecordChangeTracker();
				}
				ProgressLogger progress = new ProgressLogger(log);
				ExecutorService threadpool = Executors.newFixedThreadPool(WORKER_THREADS, new ThreadFactoryBuilder().setDaemon(false).setNameFormat("ComputeSamTags-%d").build());
    			try (SAMRecordIterator it = reader.iterator()) {
					Iterator<SAMRecord> asyncIt = transform(threadpool, Defaults.ASYNC_BUFFER_SIZE, it, tracker);
    				File tmpoutput = gridss.Defaults.OUTPUT_TO_TEMP_FILE ? FileSystemContext.getWorkingFileFor(OUTPUT, "gridss.tmp.ComputeSamTags.") : OUTPUT;
    				try (SAMFileWriter writer = writerFactory.makeSAMOrBAMWriter(header, true, tmpoutput)) {
						while (asyncIt.hasNext()) {
							SAMRecord r = asyncIt.next();
							for (String tag : REMOVE_TAGS) {
								r.setAttribute(tag, null);
							}
							writer.addAlignment(r);
							progress.record(r);
						}
    				}
    				if (tmpoutput != OUTPUT) {
    					FileHelper.move(tmpoutput, OUTPUT, true);
    				}
    			}
    			if (tracker != null) {
    				tracker.writeSummary(MODIFICATION_SUMMARY_FILE);
				}
    		}
		} catch (IOException e) {
			log.error(e);
			return -1;
		}
    	return 0;
	}
	public Iterator<SAMRecord> transform(Executor threadpool, int batchSize, Iterator<SAMRecord> it, SAMRecordChangeTracker tracker) {
		final Set<String> tags = TAGS;
		final boolean softenHardClips = SOFTEN_HARD_CLIPS;
		final boolean fixMates = FIX_MATE_INFORMATION;
		final boolean fixDuplicates = FIX_DUPLICATE_FLAG;
		final boolean fixTruncated = FIX_MISSING_HARD_CLIP;
		final boolean recalculateSupplementary = RECALCULATE_SUPPLEMENTARY;
		final int supplementaryOverlap = SUPPLEMENTARY_ALIGNMENT_OVERLAP_THRESHOLD;
		final boolean fixTerminalCigar = FIX_TERMINAL_CIGAR_INDEL;
		final ReferenceLookup reference = isReferenceRequired() ? getReference() : null;
		ParallelTransformIterator<List<SAMRecord>, List<SAMRecord>> parallelIt = new ParallelTransformIterator<>(
				SAMRecordUtil.groupedByReadName(it),
				r -> trackedTransform(
					tracker,
					r,
					tags,
					softenHardClips,
					fixMates,
					fixDuplicates,
					fixTruncated,
					recalculateSupplementary,
					supplementaryOverlap,
					fixTerminalCigar,
					reference),
				batchSize,
				threadpool);
		return new UngroupingIterator(parallelIt);
	}
	public static List<SAMRecord> trackedTransform(
			SAMRecordChangeTracker tracker,
			List<SAMRecord> records,
			Set<String> tags,
			boolean softenHardClips,
			boolean fixMates,
			boolean fixDuplicates,
			boolean fixTruncated,
			boolean recalculateSupplementary,
			int supplementaryOverlap,
			boolean fixTerminalCigar,
			ReferenceLookup reference) {
		SAMRecordChangeTracker.TrackedFragment tf = null;
		if (tracker != null) {
			tf = tracker.startTrackedChanges(records);
		}
		List<SAMRecord> after = transform(records, tags, softenHardClips, fixMates, fixDuplicates, fixTruncated, recalculateSupplementary, supplementaryOverlap, fixTerminalCigar, reference);
		if (tracker != null) {
			tracker.processTrackedChanges(tf, after);
		}
		return after;
	}
	public static List<SAMRecord> transform(
			List<SAMRecord> records,
			Set<String> tags,
			boolean softenHardClips,
			boolean fixMates,
			boolean fixDuplicates,
			boolean fixTruncated,
			boolean recalculateSupplementary,
			int supplementaryOverlap,
			boolean fixTerminalCigar,
			ReferenceLookup reference) {
		if (fixDuplicates) {
			ensureConsistentDuplicateFlag(records);
		}
		if (fixTerminalCigar) {
			for (SAMRecord r : records) {
				SAMRecordUtil.clipTerminalIndelCigars(r);
			}
		}
		if ((tags.contains(SAMTag.NM.name()) || tags.contains(SAMTag.SA.name()))) {
			for (SAMRecord r : records) {
				SAMRecordUtil.ensureNmTag(reference, r);
			}
		}
		// TODO: support secondary alignments by switching to templateBySegmentByAlignmentGroup() and handling split secondary reads
		List<SAMRecord> outRecords = new ArrayList<>(records.size());
		List<List<SAMRecord>> segments = SAMRecordUtil.templateBySegment(records);
		for (List<SAMRecord> segment : segments) {
			if (fixTruncated) {
				SAMRecordUtil.addMissingHardClipping(segment);
			}
			if (softenHardClips) {
				SAMRecordUtil.softenHardClips(segment);
			}
			// Strip out all unmapped segments
			if (anyMapped(segment)) {
				for (int i = segment.size() - 1; i >= 0; i--) {
					if (segment.get(i).getReadUnmappedFlag()) {
						segment.remove(i);
					}
				}
			}
			if (recalculateSupplementary) {
				SAMRecordUtil.reinterpretAsSplitReadAlignment(segment, supplementaryOverlap);
			} else if (tags.contains(SAMTag.SA.name())) {
				SAMRecordUtil.recalculateSupplementaryFromSA(segment);
			}
			if (Sets.intersection(tags, ImmutableSet.of(SAMTag.CC.name(), SAMTag.CP.name(), SAMTag.HI.name(), SAMTag.IH.name())).size() > 0) {
				SAMRecordUtil.calculateMultimappingTags(tags, segment);
			}
			outRecords.addAll(segment);
		}
		records = outRecords; // SAMRecords could have been removed from segments - need to reflect that in records
		if (fixMates || tags.contains(SAMTag.MC.name()) || tags.contains(SAMTag.MQ.name())) {
			SAMRecordUtil.matchReadPairPrimaryAlignments(segments);
			if (segments.size() >= 2) {
				// hack to prevent SAMRecord getters from throwing an exception
				records.stream().forEach(r -> r.setReadPairedFlag(true));
				segments.get(0).stream().forEach(r -> r.setFirstOfPairFlag(true));
				segments.get(1).stream().forEach(r -> r.setSecondOfPairFlag(true));
				SAMRecordUtil.fixMates(segments, tags.contains(SAMTag.MC.name()), tags.contains(SAMTag.MQ.name()));
			} else {
				for (SAMRecord r : records) {
					SAMRecordUtil.clearMateInformation(r, true);
				}
			}
		}
		// R2
		if (tags.contains(SAMTag.R2.name())) {
			SAMRecordUtil.calculateTagR2(segments);
		}
		// Q2
		if (tags.contains(SAMTag.Q2.name())) {
			SAMRecordUtil.calculateTagQ2(segments);
		}
		return records;
	}
	private static boolean anyMapped(List<SAMRecord> records) {
		for (SAMRecord r : records) {
			if (!r.getReadUnmappedFlag()) return true;
		}
		return false;
	}
	protected boolean isReferenceRequired() {
		return TAGS.contains(SAMTag.NM.name()) ||
			TAGS.contains(SAMTag.SA.name()) || // SA requires NM
			FIX_TERMINAL_CIGAR_INDEL;
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
