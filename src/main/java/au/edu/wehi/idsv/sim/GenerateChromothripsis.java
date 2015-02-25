package au.edu.wehi.idsv.sim;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;

import java.io.File;

import net.sf.picard.io.IoUtil;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;

import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;
import org.broadinstitute.variant.vcf.VCFHeader;

import au.edu.wehi.bioinf.appraisal.VCFStructuralVariantHeaderLines;
import au.edu.wehi.bioinf.commandline.GenerateReferenceVcf;
import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.VariantCallingParameters;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.util.BufferedReferenceSequenceFile;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.Usage;
import picard.cmdline.CommandLineProgram;

/**
 * Simulates chromothripsis through random translocation
 * @author cameron.d
 *
 */
public class GenerateChromothripsis extends CommandLineProgram {
	private Log log = Log.getInstance(GenerateChromothripsis.class);
	private static final String PROGRAM_VERSION = "0.1";

    // The following attributes define the command-line arguments
    @picard.cmdline.Usage
    public String USAGE = getStandardUsagePreamble() + "Translocation breakpoint simulator " + PROGRAM_VERSION;
    @Option(doc="Reference used for alignment", shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME)
    public File REFERENCE;
    @Option(doc="VCF variant list")
    public File VCF;
    @Option(doc="Number of bases between variants")
    public int SIZE;
    @Option(doc="Number translocations")
    public int COUNT;
    @Option(doc="Minimum of bases of unambiguous reference around breakponts")
    public int PADDING;
    @Option(doc="Seed for random number generator", optional=true)
    public int RANDOM_SEED = 1;
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
        			false);
        	// load variants
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
