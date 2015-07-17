package au.edu.wehi.idsv;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.io.IOException;

import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import au.edu.wehi.idsv.bed.BedpeWriter;
import au.edu.wehi.idsv.util.FileHelper;

import com.google.common.collect.Lists;

@CommandLineProgramProperties(
        usage = "Converts VCF breakend calls to BEDPE format. "
        		+ "All variants, including structural variations, that are not in breakend format are ignored. "
        		+ "Gridss breakpoint fields, if present, are stored in the optional columns. ",  
        usageShort = "Converts VCF breakend calls to BEDPE format."
)
public class VcfBreakendToBedpe extends picard.cmdline.CommandLineProgram {
	private Log log = Log.getInstance(VcfBreakendToBedpe.class);
	@Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="VCF containing structural variation breakend calls")
    public File INPUT;
	@Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="BEDPE output file")
    public File OUTPUT;
	@Option(doc="Include VCF filtered calls in output.")
	public boolean INCLUDE_FILTERED = false;
	@Option(shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="Reference used for alignment")
    public File REFERENCE;
	@Option(doc="Inlcude header line with column names.")
	public boolean INCLUDE_HEADER = false;
	@Option(shortName="LOW", doc="Write record at breakend with lower genomic coordinate")
	public boolean INCLUDE_LOW_BREAKEND = true;
	@Option(shortName="HIGH", doc="Write record at breakend with higher genomic coordinate")
	public boolean INCLUDE_HIGH_BREAKEND = false;
	@Override
	protected int doWork() {
		if (TMP_DIR == null || TMP_DIR.size() == 0) {
			TMP_DIR = Lists.newArrayList(new File("."));
		}
		try {
			ProcessingContext pc = new ProcessingContext(new FileSystemContext(TMP_DIR.get(0), MAX_RECORDS_IN_RAM), null, null, null, null, null, null, REFERENCE, false, false);
			writeBreakpointBedpe(pc, INPUT, OUTPUT, INCLUDE_FILTERED, INCLUDE_HEADER, INCLUDE_LOW_BREAKEND, INCLUDE_HIGH_BREAKEND);
		} catch (IOException e) {
			log.error(e);
			return -1;
		}
		return 0;
	}
	public static void writeBreakpointBedpe(ProcessingContext pc, File vcf, File bedpe, boolean includeFiltered, boolean includeHeader, boolean writeLow, boolean writeHigh) throws IOException {
		if (!writeLow && !writeHigh) {
			throw new IllegalArgumentException("No breakends to be written. At least one of {LOW, HIGH} breakends should be specified");
		}
		File working = FileSystemContext.getWorkingFileFor(bedpe);
		VCFFileReader vcfReader = null;
		CloseableIterator<VariantContext> it = null; 
		BedpeWriter writer = null;
		try {
			vcfReader = new VCFFileReader(vcf, false);
			it = vcfReader.iterator();
			writer = new BedpeWriter(pc.getDictionary(), working);
			if (includeHeader) {
				writer.writeHeader();
			}
			while (it.hasNext()) {
				IdsvVariantContext variant = IdsvVariantContext.create(pc, null, it.next());
				if (variant instanceof VariantContextDirectedBreakpoint) {
					VariantContextDirectedBreakpoint bp = (VariantContextDirectedBreakpoint)variant;
					if (bp.getBreakendSummary().isLowBreakend() && writeLow) {
						if (!bp.isFiltered() || includeFiltered) {
							writer.write(bp);
						}
					}
					if (bp.getBreakendSummary().isHighBreakend() && writeHigh) {
						if (!bp.isFiltered() || includeFiltered) {
							writer.write(bp);
						}
					}
				}
			}
			writer.close();
			FileHelper.move(working, bedpe, false);
		} finally {
			CloserUtil.close(writer);
			CloserUtil.close(it);
			CloserUtil.close(vcfReader);
		}
	}
	public static void main(String[] argv) {
        System.exit(new VcfBreakendToBedpe().instanceMain(argv));
    }
}
