package au.edu.wehi.idsv;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
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
import picard.sam.FixMateInformation;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

@CommandLineProgramProperties(
        usage = "Converts VCF breakend calls to BAM format for visualisation. "
        		+ "Each breakpoint is represented as a read pair: one for each breakend. "
        		+ "All variants, including structural variations, that are not in breakend format are ignored. "
        		+ "VCF fields are not imported into the BAM",  
        usageShort = "Converts VCF breakend calls to BAM format for visualisation."
)
public class VcfBreakendToReadPair extends picard.cmdline.CommandLineProgram {
	private Log log = Log.getInstance(VcfBreakendToReadPair.class);
	@Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="VCF containing structural variation breakend calls")
    public File INPUT;
	@Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="BAM output file")
    public File OUTPUT;
	@Option(shortName="OF", doc="BAM output file of filtered calls")
    public File OUTPUT_FILTERED;
	@Option(shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="Reference genome")
    public File REFERENCE;
	@Override
	protected int doWork() {
		if (TMP_DIR == null || TMP_DIR.size() == 0) {
			TMP_DIR = Lists.newArrayList(new File("."));
		}
		try {
			GenomicProcessingContext pc = new GenomicProcessingContext(new FileSystemContext(TMP_DIR.get(0), MAX_RECORDS_IN_RAM), REFERENCE, false, null);
			writeVisualisationBam(pc, INPUT, OUTPUT, OUTPUT_FILTERED);
		} catch (IOException e) {
			log.error(e);
			return -1;
		}
		return 0;
	}
	public void writeVisualisationBam(GenomicProcessingContext pc, File vcf, File bam, File bamFiltered) throws IOException {
		File working = FileSystemContext.getWorkingFileFor(bam);
		File workingFiltered = FileSystemContext.getWorkingFileFor(bamFiltered);
		VCFFileReader vcfReader = new VCFFileReader(vcf, false);
		CloseableIterator<VariantContext> it = vcfReader.iterator();
		SAMFileWriter writer = null;
		SAMFileWriter writerFiltered = null;
		try {
			SAMFileWriterFactory factory = pc.getSamFileWriterFactory(true);
			SAMFileHeader header = pc.getBasicSamHeader();
			writer = factory.makeBAMWriter(header, false, working);
			writerFiltered =  factory.makeBAMWriter(header, false, workingFiltered);
			while (it.hasNext()) {
				IdsvVariantContext variant = IdsvVariantContext.create(pc, null, it.next());
				if (variant instanceof VariantContextDirectedBreakpoint) {
					VariantContextDirectedBreakpoint bp = (VariantContextDirectedBreakpoint)variant;
					if (bp.isFiltered()) {
						writerFiltered.addAlignment(bp.asSamRecord(header));
					} else {
						writer.addAlignment(bp.asSamRecord(header));
					}
				}
			}
			writer.close();
			writerFiltered.close();
			// Correct mate pairing since asSAMRecord() does not factor in mate anchor cigar 
			new FixMate().fix(working, bam);
			new FixMate().fix(workingFiltered, bamFiltered);
		} finally {
			CloserUtil.close(writer);
			CloserUtil.close(writerFiltered);
			CloserUtil.close(it);
			CloserUtil.close(vcfReader);
			working.delete();
			workingFiltered.delete();
		}
	}
	private class FixMate extends FixMateInformation {
		public void fix(File file, File bam) {
			TMP_DIR = VcfBreakendToReadPair.this.TMP_DIR;
			VERBOSITY = VcfBreakendToReadPair.this.VERBOSITY;
			QUIET = VcfBreakendToReadPair.this.QUIET;
			VALIDATION_STRINGENCY = VcfBreakendToReadPair.this.VALIDATION_STRINGENCY;
			COMPRESSION_LEVEL = VcfBreakendToReadPair.this.COMPRESSION_LEVEL;
			MAX_RECORDS_IN_RAM = VcfBreakendToReadPair.this.MAX_RECORDS_IN_RAM;
			CREATE_INDEX = VcfBreakendToReadPair.this.CREATE_INDEX;
			CREATE_MD5_FILE = VcfBreakendToReadPair.this.CREATE_MD5_FILE;
			REFERENCE_SEQUENCE = VcfBreakendToReadPair.this.REFERENCE_SEQUENCE;
			GA4GH_CLIENT_SECRETS = VcfBreakendToReadPair.this.GA4GH_CLIENT_SECRETS;
			
			ASSUME_SORTED = false;
			INPUT = ImmutableList.of(file);
			OUTPUT = bam;
			SORT_ORDER = SortOrder.coordinate;
			doWork();
		}
	}
	public static void main(String[] argv) {
        System.exit(new VcfBreakendToReadPair().instanceMain(argv));
    }
}
