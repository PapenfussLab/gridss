package au.edu.wehi.idsv.sim;

import java.io.File;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.ProcessingContext;

public abstract class SimulationGenerator extends CommandLineProgram {
    @Option(doc="Reference used for alignment", shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME)
    public File REFERENCE;
    @Option(doc="Variant list")
    public File VCF;
    @Option(doc="Resultant sequence")
    public File FASTA;
    @Option(doc="Minimum of bases of unambiguous reference sequence around breakpoints", optional=true)
    public int PADDING = 1500;
    @Option(doc="Seed for random number generator", optional=true)
    public int RANDOM_SEED = 1;
    @Option(doc="Chromosome to simulate")
    public String CHR;
    @Option(doc="Include reference chromosome in output as separate contig", optional=true)
    public boolean INCLUDE_REFERENCE = false;
    protected ProcessingContext getProcessingContext() {
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
    	return pc;
    }
}