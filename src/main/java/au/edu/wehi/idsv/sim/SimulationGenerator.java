package au.edu.wehi.idsv.sim;

import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.GenomicProcessingContext;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.picard.TwoBitBufferedReferenceSequenceFile;
import gridss.cmdline.ReferenceCommandLineProgram;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.FileNotFoundException;

public abstract class SimulationGenerator extends ReferenceCommandLineProgram {
    @Argument(doc="Variant list")
    public File VCF;
    @Argument(doc="Resultant sequence", shortName="FA")
    public File FASTA;
    @Argument(doc="Minimum of bases of unambiguous reference sequence around breakpoints", optional=true)
    public int UNAMBIGUOUS_MARGIN = 1500;
    @Argument(doc="Seed for random number generator", optional=true)
    public int RANDOM_SEED = 1;
    @Argument(doc="Chromosome to simulate")
    public String CHR;
    @Argument(doc="Include reference chromosome in output as separate contig", optional=true)
    public boolean INCLUDE_REFERENCE = false;
    protected GenomicProcessingContext getProcessingContext() {
        GenomicProcessingContext pc = new GenomicProcessingContext(getFileSystemContext(), REFERENCE_SEQUENCE, getReference());
    	pc.setCommandLineProgram(this);
    	return pc;
    }
}
