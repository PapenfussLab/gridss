package au.edu.wehi.idsv;

import au.edu.wehi.idsv.bed.BedpeWriter;
import au.edu.wehi.idsv.util.FileHelper;
import com.google.common.collect.Lists;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.IOException;
import java.util.Locale;

@CommandLineProgramProperties(
		summary = "Converts VCF breakend calls to BEDPE format. "
        		+ "All variants, including structural variations, that are not in breakend format are ignored. "
        		+ "Gridss breakpoint fields, if present, are stored in the optional columns. ",  
        oneLineSummary = "Converts VCF breakend calls to BEDPE format.",
        programGroup = gridss.cmdline.programgroups.DataConversion.class
)
public class VcfBreakendToBedpe extends picard.cmdline.CommandLineProgram {
	private Log log = Log.getInstance(VcfBreakendToBedpe.class);
	@Argument(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="VCF containing structural variation breakend calls")
    public File INPUT;
	@Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="BEDPE output file containing unfiltered calls")
    public File OUTPUT;
	@Argument(shortName="OF", doc="BEDPE output file of filtered calls")
    public File OUTPUT_FILTERED;
	@Argument(shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="Reference used for alignment")
    public File REFERENCE;
	@Argument(doc="Include header line with column names.")
	public boolean INCLUDE_HEADER = false;
	@Argument(shortName="LOW", doc="Write record at breakend with lower genomic coordinate")
	public boolean INCLUDE_LOW_BREAKEND = true;
	@Argument(shortName="HIGH", doc="Write record at breakend with higher genomic coordinate")
	public boolean INCLUDE_HIGH_BREAKEND = false;
	@Override
	protected int doWork() {
		log.debug("Setting language-neutral locale");
    	java.util.Locale.setDefault(Locale.ROOT);
		if (TMP_DIR == null || TMP_DIR.size() == 0) {
			TMP_DIR = Lists.newArrayList(new File("."));
		}
		try {
			GenomicProcessingContext pc = new GenomicProcessingContext(new FileSystemContext(TMP_DIR.get(0), MAX_RECORDS_IN_RAM), REFERENCE, null);
			pc.setCommandLineProgram(this);
			writeBreakpointBedpe(pc, INPUT, OUTPUT, OUTPUT_FILTERED, INCLUDE_HEADER, INCLUDE_LOW_BREAKEND, INCLUDE_HIGH_BREAKEND);
		} catch (IOException e) {
			log.error(e);
			return -1;
		}
		return 0;
	}
	public static void writeBreakpointBedpe(GenomicProcessingContext pc, File vcf, File bedpe, File bedpeFiltered, boolean includeHeader, boolean writeLow, boolean writeHigh) throws IOException {
		if (!writeLow && !writeHigh) {
			throw new IllegalArgumentException("No breakends to be written. At least one of {LOW, HIGH} breakends should be specified");
		}
		File working = gridss.Defaults.OUTPUT_TO_TEMP_FILE ? FileSystemContext.getWorkingFileFor(bedpe) : bedpe;
		File workingFiltered = gridss.Defaults.OUTPUT_TO_TEMP_FILE ? FileSystemContext.getWorkingFileFor(bedpeFiltered) : bedpeFiltered;
		VCFFileReader vcfReader = null;
		CloseableIterator<VariantContext> it = null; 
		BedpeWriter writer = null;
		BedpeWriter writerFiltered = null;
		try {
			vcfReader = new VCFFileReader(vcf, false);
			it = vcfReader.iterator();
			writer = new BedpeWriter(pc.getDictionary(), working);
			writerFiltered = new BedpeWriter(pc.getDictionary(), workingFiltered);
			if (includeHeader) {
				writer.writeHeader(true, true);
				writerFiltered.writeHeader(true, true);
			}
			while (it.hasNext()) {
				IdsvVariantContext variant = IdsvVariantContext.create(pc.getDictionary(), null, it.next());
				if (variant instanceof VariantContextDirectedBreakpoint) {
					VariantContextDirectedBreakpoint bp = (VariantContextDirectedBreakpoint)variant;
					if (bp.getBreakendSummary().isLowBreakend() && writeLow) {
						if (bp.isFiltered()) {
							writerFiltered.write(bp);
						} else {
							writer.write(bp);
						}
					}
					if (bp.getBreakendSummary().isHighBreakend() && writeHigh) {
						if (bp.isFiltered()) {
							writerFiltered.write(bp);
						} else {
							writer.write(bp);
						}
					}
				}
			}
			writer.close();
			writerFiltered.close();
		} finally {
			CloserUtil.close(writer);
			CloserUtil.close(writerFiltered);
			CloserUtil.close(it);
			CloserUtil.close(vcfReader);
		}
		if (working != bedpe) {
			FileHelper.move(working, bedpe, false);
		}
		if (workingFiltered != bedpeFiltered) {
			FileHelper.move(workingFiltered, bedpeFiltered, false);
		}
	}
	public static void main(String[] argv) {
        System.exit(new VcfBreakendToBedpe().instanceMain(argv));
    }
}
