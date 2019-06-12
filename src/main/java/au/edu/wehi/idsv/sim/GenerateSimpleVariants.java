package au.edu.wehi.idsv.sim;

import au.edu.wehi.idsv.GenomicProcessingContext;
import au.edu.wehi.idsv.vcf.SvType;
import com.google.common.collect.Lists;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;

import java.util.List;
import java.util.Locale;

@CommandLineProgramProperties(
        summary = "Create a fasta containing structural variants of the requested types. Can simulate, insertion of random sequence, deletion, inversion, and tandem duplication.",  
        oneLineSummary = "Simple structural variant simulator",
        programGroup = gridss.cmdline.programgroups.Benchmarking.class
)
public class GenerateSimpleVariants extends SimulationGenerator {

    @Argument(doc="List of variants to insert. Valid variants are {INS, DEL, INV, DUP} for novel sequence insertion, deletion, inversion, and tandem duplication", optional=true)
    public List<SvType> TYPE = Lists.newArrayList(SvType.INS, SvType.DEL, SvType.INV, SvType.DUP);
    @Argument(doc="Variant sizes", optional=true)
	public List<Integer> SIZE = Lists.newArrayList(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 16, 20, 24, 28, 32, 48, 64, 80, 96, 112, 128, 160, 192, 224, 256, 288, 320, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536);
    @Argument(doc="Number of copies of each variant (type,size) pairing to insert. Defaults to as many copies as possible ", optional=true)
    public Integer COPIES;
    //private static Log log = Log.getInstance(GenerateSimpleVariants.class);
    @Override
	protected int doWork() {
    	try {
        	java.util.Locale.setDefault(Locale.ROOT);
        	IOUtil.assertFileIsReadable(REFERENCE);
        	GenomicProcessingContext pc = getProcessingContext();
        	SimpleVariantChromosome gen = new SimpleVariantChromosome(pc, CHR, PADDING, RANDOM_SEED);
        	gen.assemble(FASTA, VCF, INCLUDE_REFERENCE, TYPE, SIZE, COPIES == null ? Integer.MAX_VALUE : COPIES);
        } catch (Exception e) {
			e.printStackTrace();
			return 1;
		}
        return 0;
    }
	public static void main(String[] argv) {
        System.exit(new GenerateSimpleVariants().instanceMain(argv));
    }
}
