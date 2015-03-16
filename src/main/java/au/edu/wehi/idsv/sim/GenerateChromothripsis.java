package au.edu.wehi.idsv.sim;

import htsjdk.samtools.util.IOUtil;
import picard.cmdline.Option;
import au.edu.wehi.idsv.ProcessingContext;

/**
 * Simulates chromothripsis through random translocation
 * @author cameron.d
 *
 */
public class GenerateChromothripsis extends SimulationGenerator {
	private static final String PROGRAM_VERSION = "0.1";

    // The following attributes define the command-line arguments
    @picard.cmdline.Usage
    public String USAGE = getStandardUsagePreamble() + "Translocation breakpoint simulator " + PROGRAM_VERSION;
    @Option(doc="Number of genomic fragment to shatter into.", optional=true)
    public Integer FRAGMENTS;
    @Option(doc="Fragment retention rate", optional=true)
    public double RETENTION_RATE = 0.25;
    protected int doWork() {
        try {
        	IOUtil.assertFileIsReadable(REFERENCE);
        	ProcessingContext pc = getProcessingContext();
        	int referenceIndex = 0;
        	int length = pc.getReference().getSequenceDictionary().getSequence(referenceIndex).getSequenceLength();
        	int fragments = (FRAGMENTS == null || FRAGMENTS <= 0) ? (length / (2 * PADDING) / 4) : FRAGMENTS ;
        	int retainedFragments = (int)(fragments * RETENTION_RATE);
        	FragmentedChromosome fc = new FragmentedChromosome(pc, CHR, PADDING, RANDOM_SEED);
        	fc.shatter(fragments);
        	fc.assemble(FASTA, VCF, retainedFragments, INCLUDE_REFERENCE);
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
