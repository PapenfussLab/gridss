package au.edu.wehi.idsv.sim;

import htsjdk.samtools.util.IOUtil;

import java.io.File;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.ProcessingContext;

/**
 * Simulates chromothripsis through random translocation
 * @author cameron.d
 *
 */
public class GenerateChromothripsis extends CommandLineProgram {
	private static final String PROGRAM_VERSION = "0.1";

    // The following attributes define the command-line arguments
    @picard.cmdline.Usage
    public String USAGE = getStandardUsagePreamble() + "Translocation breakpoint simulator " + PROGRAM_VERSION;
    @Option(doc="Reference used for alignment", shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME)
    public File REFERENCE;
    @Option(doc="VCF variant list")
    public File VCF;
    @Option(doc="Reassembled shatter chromosome")
    public File FASTA;
    @Option(doc="Number of genomic fragment to shatter into", optional=true)
    public Integer FRAGMENTS;
    @Option(doc="Fragment retention rate", optional=true)
    public double RETENTION_RATE = 0.25;
    @Option(doc="Minimum of bases of unambiguous reference around breakponts", optional=true)
    public int PADDING = 1500;
    @Option(doc="Seed for random number generator", optional=true)
    public int RANDOM_SEED = 1;
    @Option(doc="chr to shatter")
    public String CHR;
    protected int doWork() {
        try {
        	IOUtil.assertFileIsReadable(REFERENCE);
        	ProcessingContext pc = new ProcessingContext(
        			new FileSystemContext(TMP_DIR.get(0), new File("."), MAX_RECORDS_IN_RAM),
        			null,
        			null,
        			null,
        			null,
        			null,
        			null,
        			REFERENCE,
        			false,
        			true);
        	int referenceIndex = 0;
        	int length = pc.getReference().getSequenceDictionary().getSequence(referenceIndex).getSequenceLength();
        	int fragments = (FRAGMENTS == null || FRAGMENTS <= 0) ? (length / (2 * PADDING) / 4) : FRAGMENTS ;
        	int retainedFragments = (int)(fragments * RETENTION_RATE);
        	FragmentedChromosome fc = new FragmentedChromosome(pc, CHR, PADDING, RANDOM_SEED);
        	fc.shatter(fragments);
        	fc.assemble(FASTA, VCF, retainedFragments);
        } catch (Exception e) {
			e.printStackTrace();
			return 1;
		}
        return 0;
    }
	public static void main(String[] argv) {
	    System.exit(new GenerateChromothripsis().instanceMain(argv));
	}
}
