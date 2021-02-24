package gridss;

import au.edu.wehi.idsv.GenomicProcessingContext;
import au.edu.wehi.idsv.IterativeSplitReadRealigner;
import au.edu.wehi.idsv.SplitReadRealigner;
import au.edu.wehi.idsv.StreamingSplitReadRealigner;
import au.edu.wehi.idsv.alignment.BwaStreamingAligner;
import au.edu.wehi.idsv.alignment.ExternalProcessFastqAligner;
import au.edu.wehi.idsv.alignment.ExternalProcessStreamingAligner;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import gridss.cmdline.ReferenceCommandLineProgram;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.StandardOptionDefinitions;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

@CommandLineProgramProperties(
        summary = "Uses an external aligner to identify split reads by iterative alignment of soft clipped bases. "
        		+ "Existing split read alignments are left untouched.",
        oneLineSummary = "Converts soft clipped reads to split reads",
        programGroup = gridss.cmdline.programgroups.DataCleaning.class
)
public class SoftClipsToSplitReads extends ReferenceCommandLineProgram {
	private static final Log log = Log.getInstance(SoftClipsToSplitReads.class);
	public static final List<String> BWA_COMMAND_LINE = ImmutableList.of(
			"bwa",
			"mem",
			"-K", "10000000",
			"-L", "0,0",
			"-t", "%3$d",
			"%2$s",
			"%1$s");
	public static final List<String> BOWTIE2_COMMAND_LINE = ImmutableList.of("bowtie2", "--threads", "%3$d", "--local", "--mm", "--reorder", "-x", "%2$s", "-U", "%1$s");
    @Argument(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Coordinate-sorted input file", optional=false)
    public File INPUT;
    @Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output file", optional=false)
    public File OUTPUT;
	@Argument(doc="Outputs new supplementary records and primary alignments with a new position to a separate file. " +
			"Records ordering matches initial primary alignment records INPUT." +
			"If this parameter is omitted, modified records are coordinate sorted and merged into the OUTPUT file. " +
			"This parameter is useful for off-loading the sort and merge steps to an external tool with better sorting performance than htsjdk (e.g. samtools)", optional=true)
	public File OUTPUT_UNORDERED_RECORDS = null;
    @Argument(doc="Minimum bases clipped. Generally, short read aligners are not able to uniquely align sequences shorter than 18-20 bases.", optional=true)
    public int MIN_CLIP_LENGTH = 20;
    @Argument(doc="Minimum average base quality score of clipped bases. Low quality clipped bases are indicative of sequencing errors.", optional=true)
    public float MIN_CLIP_QUAL = 5;
    @Argument(doc="Indicates whether to perform split read identification on secondary read alignments.", optional=true)
    public boolean PROCESS_SECONDARY_ALIGNMENTS = false;
    @Argument(doc="Indicates whether to perform realignment on existing chimeric alignment. If true, only the primary alignment record is retained.", optional=true)
    public boolean REALIGN_EXISTING_SPLIT_READS = false;
    @Argument(doc="Indicates whether to realign the entire read, or just the soft clipped bases.", optional=true)
    public boolean REALIGN_ENTIRE_READ = false;
    @Argument(doc="Indicates whether to adjust the primary alignment position if the total edit distance can be reduced by extending or contracting the primary alignment. "
    		+ "ComputeSamTags should be rerun to correct any changes in primary alignment position if this operation is performed.", optional=true)
    public boolean READJUST_PRIMARY_ALIGNMENT_POSITION = false;
	@Argument(doc="Write the original alignment to the OA SAM tag. Only relevant when REALIGN_ENTIRE_READ is true.", optional=true)
    public boolean WRITE_OA = true;
    @Argument(doc="Number of threads to use for realignment. Defaults to number of cores available."
			+ " Note that I/O threads are not included in this worker thread count so CPU usage can be higher than the number of worker thread.",
    		shortName="THREADS")
    public int WORKER_THREADS = Runtime.getRuntime().availableProcessors();
	@Argument(doc="Which aligner to use. GRIDSS supports in-process BWA alignment, as well as external aligners", optional=true)
	public Aligner ALIGNER = Aligner.EXTERNAL;
	@Argument(doc="Number of records to buffer when performing in-process or streaming alignment. Not applicable when performing external alignment.", optional=true)
	public int ALIGNER_BATCH_SIZE = MAX_RECORDS_IN_RAM;
    @Argument(doc="Directly pipe the input and output of the aligner instead of writing to intermediate files."
			+ " The aligner must support using \"-\" as the input filename when reading from stdin."
			+ " The sort order of the input file will not be retained.", optional=true)
	public boolean ALIGNER_STREAMING = false;
    @Argument(doc="Command line arguments to run external aligner. Aligner output should be written to stdout and the records MUST match the input fastq order."
    		+ "Java argument formatting is used with %1$s being the fastq file to align, "
    		+ "%2$s the reference genome, and %3$d the number of threads to use.", optional=true)
	public List<String> ALIGNER_COMMAND_LINE = Lists.newArrayList(BWA_COMMAND_LINE);
	@Argument(doc="Base quality score to sent to aligner if quality scores are missing.", optional=true)
	public byte FALLBACK_BASE_QUALITY = 20;
	/**
	 * Which aligner to perform the alignment with
	 */
	public enum Aligner {
    	BWAMEM,
		EXTERNAL,
	}

    @Override
	protected int doWork() {
		log.debug("Setting language-neutral locale");
    	java.util.Locale.setDefault(Locale.ROOT);
    	validateParameters();
    	GenomicProcessingContext pc = new GenomicProcessingContext(getFileSystemContext(), REFERENCE_SEQUENCE, getReference());
    	pc.setCommandLineProgram(this);
    	pc.setFilterDuplicates(IGNORE_DUPLICATES);
    	List<Closeable> toClose = new ArrayList<>();
    	SplitReadRealigner realigner;
    	try {
    		SamReaderFactory readerFactory = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE);
        	SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
        	switch (ALIGNER) {
				case BWAMEM:
					BwaStreamingAligner bwaAligner = new BwaStreamingAligner(REFERENCE_SEQUENCE, getReference().getSequenceDictionary(), WORKER_THREADS, ALIGNER_BATCH_SIZE * 150);
					realigner = new StreamingSplitReadRealigner(pc, bwaAligner, ALIGNER_BATCH_SIZE);
					toClose.add(bwaAligner);
					break;
				case EXTERNAL:
				default:
					if (ALIGNER_STREAMING) {
						ExternalProcessStreamingAligner streamingAligner = new ExternalProcessStreamingAligner(readerFactory, ALIGNER_COMMAND_LINE, REFERENCE_SEQUENCE, WORKER_THREADS, getReference().getSequenceDictionary());
						toClose.add(streamingAligner);
						realigner = new StreamingSplitReadRealigner(pc, streamingAligner, ALIGNER_BATCH_SIZE);
					} else {
						ExternalProcessFastqAligner externalAligner = new ExternalProcessFastqAligner(readerFactory, writerFactory, ALIGNER_COMMAND_LINE);
						realigner = new IterativeSplitReadRealigner(pc, externalAligner);
					}
					break;
			}
			realigner.setFallbackBaseQuality(FALLBACK_BASE_QUALITY);
			realigner.setMinSoftClipLength(MIN_CLIP_LENGTH);
			realigner.setMinSoftClipQuality(MIN_CLIP_QUAL);
			realigner.setProcessSecondaryAlignments(PROCESS_SECONDARY_ALIGNMENTS);
			realigner.setRealignExistingSplitReads(REALIGN_EXISTING_SPLIT_READS);
			realigner.setRealignEntireRecord(REALIGN_ENTIRE_READ);
			realigner.setWorkerThreads(WORKER_THREADS);
			realigner.setAdjustPrimaryAlignment(READJUST_PRIMARY_ALIGNMENT_POSITION);
			realigner.setWriteOATag(WRITE_OA);
			realigner.createSupplementaryAlignments(INPUT, OUTPUT, OUTPUT_UNORDERED_RECORDS);

			for (Closeable c : toClose) {
				c.close();
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
    	if (OUTPUT_UNORDERED_RECORDS != null) {
			IOUtil.assertFileIsWritable(OUTPUT_UNORDERED_RECORDS);
		}
	}

	public static void main(String[] argv) {
        System.exit(new SoftClipsToSplitReads().instanceMain(argv));
    }
}
