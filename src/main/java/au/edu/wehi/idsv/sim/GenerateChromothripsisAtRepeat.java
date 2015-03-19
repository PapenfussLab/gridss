package au.edu.wehi.idsv.sim;

import htsjdk.samtools.util.IOUtil;

import java.io.File;

import picard.cmdline.Option;
import au.edu.wehi.idsv.ProcessingContext;

/**
 * Simulates chromothripsis through random translocation
 * @author cameron.d
 *
 */
public class GenerateChromothripsisAtRepeat extends SimulationGenerator {
	private static final String PROGRAM_VERSION = "0.1";

    // The following attributes define the command-line arguments
    @picard.cmdline.Usage
    public String USAGE = getStandardUsagePreamble() + "Translocation breakpoint simulator - one side of each breakpoint is always in repetative sequence" + PROGRAM_VERSION;
    @Option(doc="Number of genomic fragments to assemble", optional=true)
    public Integer FRAGMENTS = 1000;
    @Option(doc="Size of each fragment", optional=true)
    public Integer FRAGMENT_SIZE = 2000;
    @Option(doc="Uncompressed RepeatMasker output file", shortName="RM")
    public File REPEATMASKER_OUTPUT;
    @Option(doc="Repeat class/family as output by repeatmasker", shortName="CF")
    public String CLASS_FAMILY = "SINE/Alu";
    protected int doWork() {
        try {
        	IOUtil.assertFileIsReadable(REFERENCE);
        	IOUtil.assertFileIsReadable(REPEATMASKER_OUTPUT);
        	ProcessingContext pc = getProcessingContext();
        	RepeatFragmentedChromosome fc = new RepeatFragmentedChromosome(pc, CHR, PADDING, RANDOM_SEED);
        	fc.assembleSingleSidedRepeats(FASTA, VCF, INCLUDE_REFERENCE,
        			FRAGMENTS, FRAGMENT_SIZE,
        			REPEATMASKER_OUTPUT, CLASS_FAMILY);
        } catch (Exception e) {
			e.printStackTrace();
			return 1;
		}
        return 0;
    }
	public static void main(String[] argv) {
	    System.exit(new GenerateChromothripsisAtRepeat().instanceMain(argv));
	}
}
