package gridss;

import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.GenomicProcessingContext;
import au.edu.wehi.idsv.StreamingSplitReadRealigner;
import au.edu.wehi.idsv.alignment.BwaStreamingAligner;
import au.edu.wehi.idsv.alignment.StreamingAligner;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import au.edu.wehi.idsv.util.FileHelper;
import au.edu.wehi.idsv.util.GroupingIterator;
import au.edu.wehi.idsv.util.UngroupingIterator;
import com.google.common.collect.Ordering;
import com.google.common.collect.Sets;
import gridss.cmdline.ReferenceCommandLineProgram;
import htsjdk.samtools.*;
import htsjdk.samtools.util.AsyncReadTaskRunner;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.checkerframework.checker.units.qual.C;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import java.util.Set;

@CommandLineProgramProperties(
		summary = "Preprocesses a readname sorted file for use in the GRIDSS assembler. " +
				" This command combines ComputeSamTags and SoftClipsToSplitReads, only exposing parameters consistent" +
				" with GRIDSS assembler input requirements. Realignment is done in-process." +
				" To better control CPU usage, ComputeSamTags uses the htsjdk non-blocking thread pool to perform it's work." +
				" Use -Dsamjdk.async_io_read_threads to control thread count.",
        oneLineSummary = "Preprocesses a readname sorted file for use in the GRIDSS assembler.",
        programGroup = picard.cmdline.programgroups.ReadDataManipulationProgramGroup.class
)
public class PreprocessForBreakendAssembly extends ReferenceCommandLineProgram {
	private static final Log log = Log.getInstance(PreprocessForBreakendAssembly.class);
	@Argument(doc="Minimum bases clipped. Generally, short read aligners are not able to uniquely align sequences shorter than 18-20 bases.", optional=true)
	public int MIN_CLIP_LENGTH = new SoftClipsToSplitReads().MIN_CLIP_LENGTH;
	@Argument(doc="Minimum average base quality score of clipped bases. Low quality clipped bases are indicative of sequencing errors.", optional=true)
	public float MIN_CLIP_QUAL = new SoftClipsToSplitReads().MIN_CLIP_QUAL;
	@Argument(doc="Which in-process aligner to use.", optional=true)
	public SoftClipsToSplitReads.Aligner ALIGNER = SoftClipsToSplitReads.Aligner.BWAMEM;
	@Argument(doc="Number of records to buffer when performing in-process or streaming alignment. Not applicable when performing external alignment.", optional=true)
	public int ALIGNER_BATCH_SIZE = 100000;
	@Argument(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input BAM file grouped by read name.")
	public File INPUT;
	@Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Unsorted BAM file with tags corrected and split reads identified.")
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
	// ComputeSamTags arguments that we require for GRIDSS
	private boolean FIX_DUPLICATE_FLAG = true;
	private boolean FIX_MATE_INFORMATION = true;
	private boolean FIX_SA = true;
	private boolean FIX_MISSING_HARD_CLIP = true;
	private boolean RECALCULATE_SA_SUPPLEMENTARY = true;
	private boolean SOFTEN_HARD_CLIPS = true;
	@Argument(doc="Number of threads to use for realignment. Defaults to number of cores available."
			+ " Note that I/O threads are not included in this worker thread count so CPU usage can be higher than the number of worker thread.",
			shortName="THREADS")
	public int WORKER_THREADS = Runtime.getRuntime().availableProcessors();
	@Argument(doc="Base quality score to sent to aligner if quality scores are missing.", optional=true)
	public byte FALLBACK_BASE_QUALITY = new SoftClipsToSplitReads().FALLBACK_BASE_QUALITY;
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
		return super.customCommandLineValidation();
	}

	@Override
	protected int doWork() {
		log.debug("Setting language-neutral locale");
		java.util.Locale.setDefault(Locale.ROOT);
		// validate parameters
		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsWritable(OUTPUT);
		IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);

		GenomicProcessingContext pc = new GenomicProcessingContext(getFileSystemContext(), REFERENCE_SEQUENCE, getReference());
		StreamingAligner sa;
		StreamingSplitReadRealigner realigner;
		switch (ALIGNER) {
			case BWAMEM:
				// Our bwa interface uses a basepair-based buffer size, not a record-based buffer size
				// If we choose a number too big, we'll only ever actually invoke bwa when we flush
				// which will halve our throughput.
				// 25bp per read ensures we're unlikely to be forced to flush.
				int bwaBufferSizeInBases = ALIGNER_BATCH_SIZE * 25;
				sa = new BwaStreamingAligner(REFERENCE_SEQUENCE, getReference().getSequenceDictionary(), WORKER_THREADS, bwaBufferSizeInBases);
				break;
			case EXTERNAL:
			default:
				throw new IllegalArgumentException("Aligner not supported by PreprocessForBreakendAssembly");
		}
		realigner = new StreamingSplitReadRealigner(pc, sa, ALIGNER_BATCH_SIZE);
		realigner.setFallbackBaseQuality(FALLBACK_BASE_QUALITY);
		realigner.setMinSoftClipLength(MIN_CLIP_LENGTH);
		realigner.setMinSoftClipQuality(MIN_CLIP_QUAL);
		realigner.setWorkerThreads(WORKER_THREADS);
		realigner.setProcessSecondaryAlignments(false);
		realigner.setRealignExistingSplitReads(false);
		realigner.setAdjustPrimaryAlignment(false);
		realigner.setRealignEntireRecord(false);
		realigner.setRealignExistingSplitReads(false);
		ComputeSamTags tags = new ComputeSamTags();
		tags.setReference(this.getReference());
		tags.WORKER_THREADS = -1; // don't use ComputeSamTags workers - we're handling this ourselves
		tags.TAGS = TAGS;
		tags.FIX_DUPLICATE_FLAG = FIX_DUPLICATE_FLAG;
		tags.FIX_MATE_INFORMATION = FIX_MATE_INFORMATION;
		tags.FIX_SA = FIX_SA;
		tags.FIX_MISSING_HARD_CLIP = FIX_MISSING_HARD_CLIP;
		tags.RECALCULATE_SA_SUPPLEMENTARY = RECALCULATE_SA_SUPPLEMENTARY;
		tags.SOFTEN_HARD_CLIPS = SOFTEN_HARD_CLIPS;
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
						// We can reuse the non-block task thread pool since the transforms aren't blocking operations
						Iterator<SAMRecord> tagFixedIt = tags.transform(AsyncReadTaskRunner.getNonBlockingThreadpool(), Defaults.ASYNC_BUFFER_SIZE, asyncIn);
						realigner.process(tagFixedIt, writer, writer);
					}
					if (tmpOutput != OUTPUT) {
						FileHelper.move(tmpOutput, OUTPUT, true);
					}
				}
			}
			sa.close();
		} catch (IOException e) {
			log.error(e);
			return -1;
		}
		return 0;
	}
}
