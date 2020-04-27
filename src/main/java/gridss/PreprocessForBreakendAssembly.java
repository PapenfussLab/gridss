package gridss;

import au.edu.wehi.idsv.*;
import au.edu.wehi.idsv.alignment.BwaStreamingAligner;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import au.edu.wehi.idsv.util.FileHelper;
import au.edu.wehi.idsv.util.GroupingIterator;
import au.edu.wehi.idsv.util.UngroupingIterator;
import com.google.common.collect.Ordering;
import com.google.common.collect.Sets;
import gridss.cmdline.ReferenceCommandLineProgram;
import htsjdk.samtools.*;
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
import java.util.List;
import java.util.Locale;
import java.util.Set;
import java.util.stream.Stream;

@CommandLineProgramProperties(
		summary = "Preprocesses a readname sorted file for use in the GRIDSS assembler. " +
				" This command combines ComputeSamTags and SoftClipsToSplitReads, only exposing parameters consistent" +
				" with GRIDSS assembler input requirements. Realignment is done in-process.",
        oneLineSummary = "Preprocesses a readname sorted file for use in the GRIDSS assembler.",
        programGroup = picard.cmdline.programgroups.ReadDataManipulationProgramGroup.class
)
public class PreprocessForBreakendAssembly extends ReferenceCommandLineProgram {
	private static final Log log = Log.getInstance(PreprocessForBreakendAssembly.class);
	public int MIN_CLIP_LENGTH = 15;
	@Argument(doc="Minimum average base quality score of clipped bases. Low quality clipped bases are indicative of sequencing errors.", optional=true)
	public float MIN_CLIP_QUAL = 5;
	@Argument(doc="Which aligner to use. GRIDSS supports in-process BWA alignment, as well as external aligners", optional=true)
	public SoftClipsToSplitReads.Aligner ALIGNER = SoftClipsToSplitReads.Aligner.BWAMEM;
	@Argument(doc="Number of records to buffer when performing in-process or streaming alignment. Not applicable when performing external alignment.", optional=true)
	public int ALIGNER_BATCH_SIZE = 10000;
	@Argument(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input BAM file grouped by read name.")
	public File INPUT;
	@Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Annotated BAM file.")
	public File OUTPUT;
	@Argument(shortName=StandardOptionDefinitions.ASSUME_SORTED_SHORT_NAME, doc="Assume that all records with the same read name are consecutive. "
			+ "Incorrect tags will be written if this is not the case.", optional=true)
	public boolean ASSUME_SORTED = false;
	@Argument(shortName="T", doc="Tags to calculate")
	public Set<String> TAGS = Sets.newHashSet(
			SAMTag.NM.name(),
			SAMTag.SA.name(),
			SAMTag.R2.name(),
			//SAMTag.Q2.name(), // dropping Q2 to improve runtime performance and file size
			SAMTag.MC.name(),
			SAMTag.MQ.name());

	@Argument(doc="Number of threads to use for realignment. Defaults to number of cores available."
			+ " Note that I/O threads are not included in this worker thread count so CPU usage can be higher than the number of worker thread.",
			shortName="THREADS")
	public int WORKER_THREADS = Runtime.getRuntime().availableProcessors();

	public static void main(String[] argv) {
        System.exit(new PreprocessForBreakendAssembly().instanceMain(argv));
    }

	@Override
	protected String[] customCommandLineValidation() {
		if (ALIGNER == ALIGNER.EXTERNAL) {
			return new String[]{"Cannot use external aligner for PreprocessForBreakendAssembly. Run ComputeSamTags and SoftClipsToSplitReads separately instead."};
		}
		if (WORKER_THREADS < 1) {
			return new String[]{"WORKER_THREADS must be at least 1."};
		}
		for (SAMTag tag : new SAMTag[] { SAMTag.R2, SAMTag.NM, SAMTag.SA, SAMTag.MC, SAMTag.MQ}) {
			if (!TAGS.contains(tag.name())) {
				return new String[]{"TAGS must contain " + tag.name()};
			}
		}
		return new String[]{};
	}

	@Override
	protected int doWork() {
		log.debug("Setting language-neutral locale");
		java.util.Locale.setDefault(Locale.ROOT);
		// validate parameters
		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsWritable(OUTPUT);
		// set up
		ComputeSamTags cst = new ComputeSamTags();
		copyInputs(cst);
		cst.TAGS = TAGS;
		cst.INPUT = INPUT;
		cst.OUTPUT = OUTPUT;
		cst.ASSUME_SORTED = ASSUME_SORTED;
		cst.SOFTEN_HARD_CLIPS = true;
		cst.TAGS = TAGS;
		cst.FIX_DUPLICATE_FLAG = true;
		cst.FIX_MATE_INFORMATION = true;
		cst.FIX_SA = true;
		cst.FIX_MISSING_HARD_CLIP = true;
		cst.RECALCULATE_SA_SUPPLEMENTARY = true;
		cst.WORKING_DIR = WORKING_DIR;
		cst.validateParameters();

		GenomicProcessingContext pc = new GenomicProcessingContext(getFileSystemContext(), REFERENCE_SEQUENCE, getReference());
		StreamingSplitReadRealigner realigner;
		switch (ALIGNER) {
			case BWAMEM:
				BwaStreamingAligner bwaAligner = new BwaStreamingAligner(REFERENCE_SEQUENCE, getReference().getSequenceDictionary(), WORKER_THREADS, ALIGNER_BATCH_SIZE * 150);
				realigner = new StreamingSplitReadRealigner(pc, bwaAligner, ALIGNER_BATCH_SIZE);
				break;
			case EXTERNAL:
			default:
				throw new IllegalArgumentException("Aligner not supported by PreprocessForBreakendAssembly");
		}
		realigner.setMinSoftClipLength(MIN_CLIP_LENGTH);
		realigner.setMinSoftClipQuality(MIN_CLIP_QUAL);
		realigner.setWorkerThreads(WORKER_THREADS);
		realigner.setProcessSecondaryAlignments(false);
		realigner.setRealignExistingSplitReads(false);
		realigner.setAdjustPrimaryAlignment(false);
		realigner.setRealignEntireRecord(false);
		realigner.setRealignAnchoringBases(false);
		realigner.setRealignExistingSplitReads(false);
		ProgressLogger progress = new ProgressLogger(log);
		String threadPrefix = INPUT.getName() + "-";



		try {
			SamReaderFactory readerFactory = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE);
			SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
			try (SamReader reader = readerFactory.open(INPUT)) {
				SAMFileHeader header = reader.getFileHeader();
				if (!ASSUME_SORTED) {
					if (header.getSortOrder() != SAMFileHeader.SortOrder.queryname) {
						log.error("INPUT is not sorted by queryname. "
								+ "ComputeSamTags requires that reads with the same name be sorted together. "
								+ "If the input file satisfies this constraint (the output from many aligners do),"
								+ " this check can be disabled with the ASSUME_SORTED option.");
						return -1;
					}
				}
				// strip header because our output is unordered
				header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
				try (SAMRecordIterator it = reader.iterator()) {
					File tmpOutput = gridss.Defaults.OUTPUT_TO_TEMP_FILE ? FileSystemContext.getWorkingFileFor(OUTPUT, "gridss.tmp.PreprocessForBReakendAssembly.") : OUTPUT;
					try (SAMFileWriter writer = writerFactory.makeSAMOrBAMWriter(header, true, tmpOutput)) {
						CloseableIterator<SAMRecord> asyncIn = new AsyncBufferedIterator<>(it, threadPrefix + "raw");
						Iterator<SAMRecord> perRecord = cst.computePerRecordFields(asyncIn, cst.getReference(), cst.TAGS, cst.SOFTEN_HARD_CLIPS, cst.FIX_MATE_INFORMATION, cst.FIX_DUPLICATE_FLAG, cst.FIX_SA, cst.FIX_MISSING_HARD_CLIP, cst.RECALCULATE_SA_SUPPLEMENTARY);
						CloseableIterator<SAMRecord> asyncPerRecord = new AsyncBufferedIterator<>(perRecord, threadPrefix + "nm");
						Iterator<List<SAMRecord>> groupByFragment = new GroupingIterator(asyncPerRecord, Ordering.natural().onResultOf((SAMRecord r) -> r.getReadName()));
						Iterator<List<SAMRecord>> perFragment = ComputeSamTags.computePerFragmentFields(groupByFragment, cst.getReference(), cst.TAGS, cst.SOFTEN_HARD_CLIPS, cst.FIX_MATE_INFORMATION, cst.FIX_DUPLICATE_FLAG, cst.FIX_SA, cst.FIX_MISSING_HARD_CLIP, cst.RECALCULATE_SA_SUPPLEMENTARY);
						CloseableIterator<List<SAMRecord>> asyncPerFragment = new AsyncBufferedIterator<>(perFragment, threadPrefix + "tags");
						Iterator<SAMRecord> computedPerRead = new UngroupingIterator<>(asyncPerFragment);
						realigner.process(computedPerRead, writer, writer);
						writer.close();
					}
					if (tmpOutput != OUTPUT) {
						FileHelper.move(tmpOutput, OUTPUT, true);
					}
				}
			}
		} catch (IOException e) {
			log.error(e);
			return -1;
		}
		return 0;
	}
}
