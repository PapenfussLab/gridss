package gridss;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Locale;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.GenomicProcessingContext;
import au.edu.wehi.idsv.SplitReadRealigner;
import au.edu.wehi.idsv.alignment.ExternalProcessFastqAligner;
import au.edu.wehi.idsv.alignment.FastqAligner;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(
        usage = "Uses an external aligner to identify split reads by iterative alignment of soft clipped bases.",
        usageShort = "Converts soft clipped reads to split reads"
)
public class SoftClipsToSplitReads extends CommandLineProgram {
	public static final List<String> BWA_COMMAND_LINE = ImmutableList.of("bwa", "mem", "-t", "%3$d", "-M", "%2$s", "%1$s");
	public static final List<String> BOWTIE2_COMMAND_LINE = ImmutableList.of("bowtie2", "--threads", "%3$d", "--local", "--mm", "--reorder", "-x", "%2$s", "-U", "%1$s");
	private static final Log log = Log.getInstance(SoftClipsToSplitReads.class);
    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input file", optional=false)
    public File INPUT;
    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output file", optional=false)
    public File OUTPUT;
    @Option(shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="Reference genome", optional=false)
    public File REFERENCE;
    @Option(doc="Minimum bases clipped. Generally, short read aligners are not able to align sequences shorter than 18-20 bases.", optional=true)
    public int MIN_CLIP_LENGTH = 15;
    @Option(doc="Command line arguments to run external aligner. Aligner output should be written to stdout and the records MUST match the input fastq order."
    		+ "Java argument formatting is used with %1$s being the fastq file to align, "
    		+ "%2$s the reference genome, and %3$d the number of alignment threads to use.", optional=true)
    public List<String> ALIGNER_COMMAND_LINE = Lists.newArrayList(BWA_COMMAND_LINE);
    private FastqAligner createAligner() {
    	SamReaderFactory readerFactory = SamReaderFactory.make();
    	SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
    	FastqAligner aligner = new ExternalProcessFastqAligner(readerFactory, writerFactory, ALIGNER_COMMAND_LINE);
    	return aligner;
    }
    @Override
	protected int doWork() {
		log.debug("Setting language-neutral locale");
    	java.util.Locale.setDefault(Locale.ROOT);
    	validateParameters();
    	
    	FastqAligner aligner = createAligner();
    	GenomicProcessingContext pc = new GenomicProcessingContext(new FileSystemContext(TMP_DIR.get(0), new File("."), MAX_RECORDS_IN_RAM), REFERENCE, false, null);
    	SplitReadRealigner realigner = new SplitReadRealigner(pc, aligner);
    	try {
    		realigner.createSupplementaryAlignments(INPUT, OUTPUT);
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
        System.exit(new SoftClipsToSplitReads().instanceMain(argv));
    }
}
