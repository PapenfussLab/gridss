package au.edu.wehi.idsv.sim;

import htsjdk.samtools.util.IOUtil;

import java.util.List;

import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.vcf.SvType;

import com.google.common.collect.Lists;

@CommandLineProgramProperties(
        usage = "Create a fasta containing structural variants of the requested types. Can simulate, insertion of random sequence, deletion, inversion, and tandem duplication.",  
        usageShort = "Simple structural variant simulator"
)
public class GenerateSimpleVariants extends SimulationGenerator {

    @Option(doc="List of variants to insert. Valid variants are {INS, DEL, INV, DUP} for novel sequence insertion, deletion, inversion, and tandem duplication", optional=true)
    public List<SvType> TYPE = Lists.newArrayList(SvType.INS, SvType.DEL, SvType.INV, SvType.DUP);
    @Option(doc="Variant sizes", optional=true)
	public List<Integer> SIZE = Lists.newArrayList(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 16, 20, 24, 28, 32, 48, 64, 80, 96, 112, 128, 160, 192, 224, 256, 288, 320, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536);
    @Option(doc="Number of copies of each variant (type,size) pairing to insert. Defaults to as many copies as possible ", optional=true)
    public Integer COPIES;
    //private static Log log = Log.getInstance(GenerateSimpleVariants.class);
    @Override
	protected int doWork() {
    	try {
        	IOUtil.assertFileIsReadable(REFERENCE);
        	ProcessingContext pc = getProcessingContext();
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
