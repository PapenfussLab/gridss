package au.edu.wehi.idsv.sim;

import au.edu.wehi.idsv.GenomicProcessingContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;

import java.io.File;
import java.util.Locale;

/**
 * Simulates chromothripsis through random translocation
 * @author Daniel Cameron
 *
 */
@CommandLineProgramProperties(
        summary = "Shatters a chromosome into fragments of the given size, and randomly reassembles a subset of them.",  
        oneLineSummary = "Simulates chromothripsis through random translocation",
        programGroup = gridss.cmdline.programgroups.Benchmarking.class
)
public class GenerateChromothripsis extends SimulationGenerator {
    @Argument(doc="Number of genomic fragments to assemble", optional=true)
    public Integer FRAGMENTS = 1000;
    @Argument(doc="Minimum size of each fragment", optional=true)
    public Integer MIN_FRAGMENT_SIZE = 2000;
	@Argument(doc="Maximum size of each fragment", optional=true)
	public Integer MAX_FRAGMENT_SIZE = 2000;
	@Argument(doc="Length of telomeric sequence to retain at the start/end of the rearranged chromosome", optional=true)
	public Integer TELOMERE_LENGTH = 50000;
    @Argument(doc="Uncompressed RepeatMasker output file. If a file is specified, one side of each fragment will be of the specified repeat", shortName="RM", optional=true)
    public File REPEATMASKER_OUTPUT = null;
    @Argument(doc="Repeat class/family as output by repeatmasker", shortName="CF", optional=true)
    public String CLASS_FAMILY = "SINE/Alu";
    protected int doWork() {
        try {
        	java.util.Locale.setDefault(Locale.ROOT);
        	GenomicProcessingContext pc = getProcessingContext();
        	FragmentedChromosome fc;
        	if (REPEATMASKER_OUTPUT == null) {
        		fc = new FragmentedChromosome(pc, CHR, TELOMERE_LENGTH, UNAMBIGUOUS_MARGIN, MIN_FRAGMENT_SIZE, MAX_FRAGMENT_SIZE, RANDOM_SEED);
        	} else {
        		fc = new RepeatFragmentedChromosome(pc, CHR, TELOMERE_LENGTH, UNAMBIGUOUS_MARGIN, MIN_FRAGMENT_SIZE, MAX_FRAGMENT_SIZE, REPEATMASKER_OUTPUT, CLASS_FAMILY, RANDOM_SEED);
        	}
        	fc.assemble(FASTA, VCF, FRAGMENTS, INCLUDE_REFERENCE);
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
