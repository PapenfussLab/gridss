package gridss;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Locale;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

import au.edu.wehi.idsv.GenomicProcessingContext;
import au.edu.wehi.idsv.SplitReadRealigner;
import au.edu.wehi.idsv.alignment.ExternalProcessFastqAligner;
import au.edu.wehi.idsv.alignment.ExternalProcessStreamingAligner;
import gridss.cmdline.ReferenceCommandLineProgram;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(
        summary = "Uses an external aligner to identify split reads by iterative alignment of soft clipped bases. "
        		+ "Existing split read alignments are left untouched.",
        oneLineSummary = "Converts soft clipped reads to split reads",
        programGroup = gridss.cmdline.programgroups.DataCleaning.class
)
public class SoftClipsToSplitReads extends ReferenceCommandLineProgram {
	private static final Log log = Log.getInstance(SoftClipsToSplitReads.class);
	public static final List<String> BWA_COMMAND_LINE = ImmutableList.of("bwa", "mem", "-t", "%3$d", "%2$s", "%1$s");
	public static final List<String> BOWTIE2_COMMAND_LINE = ImmutableList.of("bowtie2", "--threads", "%3$d", "--local", "--mm", "--reorder", "-x", "%2$s", "-U", "%1$s");
    @Argument(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input file", optional=false)
    public File INPUT;
    @Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output file", optional=false)
    public File OUTPUT;
    @Argument(doc="Minimum bases clipped. Generally, short read aligners are not able to uniquely align sequences shorter than 18-20 bases.", optional=true)
    public int MIN_CLIP_LENGTH = 15;
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
    public boolean READJUST_PRIMARY_ALIGNMENT_POSITON = false;
    @Argument(doc="Number of threads to use for realignment. Defaults to number of cores available."
			+ " Note that I/O threads are not included in this worker thread count so CPU usage can be higher than the number of worker thread.",
    		shortName="THREADS")
    public int WORKER_THREADS = Runtime.getRuntime().availableProcessors();
    @Argument(doc="Directly pipe the input and output of the aligner instead of writing to intermediate files."
			+ " The aligner must support using \"-\" as the input filename when reading from stdin."
			+ " The sort order of the input file will not be retained.", optional=true)
	public boolean ALIGNER_STREAMING = false;
    @Argument(doc="Command line arguments to run external aligner. Aligner output should be written to stdout and the records MUST match the input fastq order."
    		+ "Java argument formatting is used with %1$s being the fastq file to align, "
    		+ "%2$s the reference genome, and %3$d the number of threads to use.", optional=true)
    public List<String> ALIGNER_COMMAND_LINE = Lists.newArrayList(BWA_COMMAND_LINE);
    @Override
	protected int doWork() {
		log.debug("Setting language-neutral locale");
    	java.util.Locale.setDefault(Locale.ROOT);
    	validateParameters();
    	GenomicProcessingContext pc = new GenomicProcessingContext(getFileSystemContext(), REFERENCE_SEQUENCE, getReference());
    	pc.setCommandLineProgram(this);
    	pc.setFilterDuplicates(IGNORE_DUPLICATES);
    	SplitReadRealigner realigner = new SplitReadRealigner(pc);
    	realigner.setMinSoftClipLength(MIN_CLIP_LENGTH);
    	realigner.setMinSoftClipQuality(MIN_CLIP_QUAL);
    	realigner.setProcessSecondaryAlignments(PROCESS_SECONDARY_ALIGNMENTS);
    	realigner.setRealignExistingSplitReads(REALIGN_EXISTING_SPLIT_READS);
    	realigner.setRealignEntireRecord(REALIGN_ENTIRE_READ);
    	realigner.setWorkerThreads(WORKER_THREADS);
    	realigner.setReference(getReference());
    	realigner.setAdjustPrimaryAlignment(READJUST_PRIMARY_ALIGNMENT_POSITON);
    	try {
    		SamReaderFactory readerFactory = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE);
        	SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
        	
        	if (ALIGNER_STREAMING) {
        		ExternalProcessStreamingAligner aligner = new ExternalProcessStreamingAligner(readerFactory, ALIGNER_COMMAND_LINE, REFERENCE_SEQUENCE, WORKER_THREADS);
        		realigner.createSupplementaryAlignments(aligner, INPUT, OUTPUT, MAX_RECORDS_IN_RAM);
        	} else {
        		ExternalProcessFastqAligner aligner = new ExternalProcessFastqAligner(readerFactory, writerFactory, ALIGNER_COMMAND_LINE);
        		realigner.createSupplementaryAlignments(aligner, INPUT, OUTPUT);
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
	@Override
	protected String[] customCommandLineValidation() {
		return super.customCommandLineValidation();
	}
	public static void main(String[] argv) {
        System.exit(new SoftClipsToSplitReads().instanceMain(argv));
    }
}
