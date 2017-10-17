package au.edu.wehi.idsv.sim;

import java.io.File;

import org.broadinstitute.barclay.argparser.Argument;

import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.GenomicProcessingContext;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

public abstract class SimulationGenerator extends CommandLineProgram {
    @Argument(doc="Reference used for alignment", shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME)
    public File REFERENCE;
    @Argument(doc="Variant list")
    public File VCF;
    @Argument(doc="Resultant sequence")
    public File FASTA;
    @Argument(doc="Minimum of bases of unambiguous reference sequence around breakpoints", optional=true)
    public int PADDING = 1500;
    @Argument(doc="Seed for random number generator", optional=true)
    public int RANDOM_SEED = 1;
    @Argument(doc="Chromosome to simulate")
    public String CHR;
    @Argument(doc="Include reference chromosome in output as separate contig", optional=true)
    public boolean INCLUDE_REFERENCE = false;
    protected GenomicProcessingContext getProcessingContext() {
    	GenomicProcessingContext pc = new GenomicProcessingContext(new FileSystemContext(TMP_DIR.get(0), new File("."), MAX_RECORDS_IN_RAM), REFERENCE, null);
    	pc.setCommandLineProgram(this);
    	return pc;
    }
}
