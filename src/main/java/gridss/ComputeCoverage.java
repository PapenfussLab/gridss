package gridss;

import java.io.File;
import java.util.Iterator;

import com.google.common.collect.Iterators;

import au.edu.wehi.idsv.CoverageCalculationMethod;
import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.GcBiasAdjuster;
import au.edu.wehi.idsv.GenomicProcessingContext;
import au.edu.wehi.idsv.IdsvVariantContext;
import au.edu.wehi.idsv.IntervalCoverageAccumulator;
import au.edu.wehi.idsv.ReadGcSummary;
import au.edu.wehi.idsv.VariantContextDirectedEvidence;
import gridss.cmdline.GcSinglePassSamProgram;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(
        usage = "Computes reference genome coverage for a given BAM", usageShort = "Compute coverage"
)
public class ComputeCoverage extends GcSinglePassSamProgram {
	@Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input BAM file grouped by read name.")
    public File INPUT;
	@Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Coverage BED")
    public File OUTPUT;
	@Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="GC-adjusted coverage BED", optional=true)
    public File GC_OUTPUT;
	@Option(shortName="B", doc="Bin size used to output coverage", optional=true)
    public int BIN_SIZE = 100;
	@Option(shortName="V", doc="GRIDSS VCF containing breakpoint to split bins at", optional=true)
	public File VCF;
	@Option(shortName="GC", doc="GC adjustment file. This file must contain two tab-seperated columns without any header lines."
			+ " The first column must contain the GC percentage for adjustment, the second the adjustment multipler."
			+ " The first column should contain the integer values 0-100 inclusive, the second a floating point adjustment multiplier."
			+ " ", optional=true)
	public File GC_ADJUSTMENT;
	@Option(shortName="CM", doc="Approach used to calculate coverage. Valid values are READ, which calculates the"
			+ " actual aligned sequence coverage, and FRAGMENT which calculated physical coverage based on the"
			+ " alignment of read pairs.", optional=true)
	public CoverageCalculationMethod COVERAGE_METHOD = CoverageCalculationMethod.READ;
	
	private IntervalCoverageAccumulator ica_gc;
	private IntervalCoverageAccumulator ica_raw;
	private GcBiasAdjuster gcAdjust;
	@Override
	protected void setup(SAMFileHeader header, File samFile) {
		ica_gc = initIntervalCoverageAccumulator();
		ica_raw = initIntervalCoverageAccumulator();
	}
	private IntervalCoverageAccumulator initIntervalCoverageAccumulator() {
		SAMSequenceDictionary dictionary = getReference().getSequenceDictionary();
		if (VCF == null) {
			return new IntervalCoverageAccumulator(COVERAGE_METHOD, dictionary, BIN_SIZE, null);
		} else {
			try (VCFFileReader vcfReader = new VCFFileReader(VCF, false)) {
				try (CloseableIterator<VariantContext> it = vcfReader.iterator()) {
					GenomicProcessingContext pc = new GenomicProcessingContext(new FileSystemContext(TMP_DIR.get(0), MAX_RECORDS_IN_RAM), REFERENCE_SEQUENCE, getReference());
					Iterator<IdsvVariantContext> idsvIt = Iterators.transform(it, variant -> IdsvVariantContext.create(pc, null, variant));
					Iterator<VariantContextDirectedEvidence> bpit = Iterators.filter(idsvIt, VariantContextDirectedEvidence.class);
					return new IntervalCoverageAccumulator(COVERAGE_METHOD, dictionary, BIN_SIZE, bpit);
				}
			}
		}
	}
	@Override
	protected void acceptRead(SAMRecord record, ReferenceSequence refSeq) {
		ReadGcSummary gc = new ReadGcSummary(record, refSeq, UNPAIRED_FRAGMENT_SIZE, getReadPairConcordanceCalculator());
		ica_gc.add(record, gc, gcAdjust.adjustmentMultiplier((int)gc.gcPercentage));
		ica_raw.add(record, gc, 1.0);
	}
	@Override
	protected void finish() {
		// Write BED files
		ica_raw.writeToBed(OUTPUT);
		ica_gc.writeToBed(GC_OUTPUT);
	}
}
