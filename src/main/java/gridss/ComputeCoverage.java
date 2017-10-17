package gridss;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;

import com.google.common.collect.Iterators;

import au.edu.wehi.idsv.CoverageCalculationMethod;
import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.GcBiasAdjuster;
import au.edu.wehi.idsv.GenomicProcessingContext;
import au.edu.wehi.idsv.IdsvVariantContext;
import au.edu.wehi.idsv.IntervalCoverageAccumulator;
import au.edu.wehi.idsv.PrecomputedGcBiasAdjuster;
import au.edu.wehi.idsv.ReadGcSummary;
import au.edu.wehi.idsv.VariantContextDirectedEvidence;
import gridss.cmdline.GcSinglePassSamProgram;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

@CommandLineProgramProperties(
		summary = "Computes reference genome coverage for a given BAM",
		oneLineSummary = "Computes coverage",
		programGroup = picard.cmdline.programgroups.SamOrBam.class
)
public class ComputeCoverage extends GcSinglePassSamProgram {
	private static final Log log = Log.getInstance(ComputeCoverage.class);
	@Argument(shortName="GCO", doc="GC-adjusted coverage BED", optional=true)
    public File OUTPUT_GC;
	@Argument(shortName="B", doc="Bin size used to output coverage", optional=true)
    public int BIN_SIZE = 1000;
	@Argument(shortName="V", doc="GRIDSS VCF containing breakpoint to split bins at", optional=true)
	public File VCF;
	@Argument(shortName="GC", doc="GC adjustment file. This file must contain two tab-seperated columns without any header lines."
			+ " The first column must contain the GC percentage for adjustment, the second the adjustment multipler."
			+ " The first column should contain the integer values 0-100 inclusive, the second a floating point adjustment multiplier."
			+ " ", optional=true)
	public File GC_ADJUSTMENT;
	@Argument(shortName="CM", doc="Approach used to calculate coverage. Valid values are READ, which calculates the"
			+ " actual aligned sequence coverage, and FRAGMENT which calculated physical coverage based on the"
			+ " alignment of read pairs.", optional=true)
	public CoverageCalculationMethod COVERAGE_METHOD = CoverageCalculationMethod.READ;
	
	private IntervalCoverageAccumulator ica_gc;
	private IntervalCoverageAccumulator ica_raw;
	private GcBiasAdjuster gcAdjust;
	@Override
	protected String[] customCommandLineValidation() {
		if (OUTPUT_GC != null) {
			if (GC_ADJUSTMENT == null) {
				return new String[] { "GC_ADJUSTMENT file is required if GC_OUTPUT specified" };
			}
		}
		return super.customCommandLineValidation();
	}
	@Override
	protected void setup(SAMFileHeader header, File samFile) {
		ica_gc = null;
		if (GC_ADJUSTMENT != null) {
			try {
				gcAdjust = new PrecomputedGcBiasAdjuster(GC_ADJUSTMENT);
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
			ica_gc = initIntervalCoverageAccumulator();
		}
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
		if (ica_gc != null) {
			ica_gc.add(record, gc, gcAdjust.adjustmentMultiplier((int)gc.gcPercentage));
		}
		ica_raw.add(record, gc, 1.0);
	}
	@Override
	protected void finish() {
		// Write BED files
		try {
			ica_raw.writeToBed(OUTPUT);
		} catch (IOException e) {
			String msg = String.format("Unable to write to %s", OUTPUT);
			log.error(e, msg);
			throw new RuntimeException(msg, e);
		}
		if (ica_gc != null) {
			try {
				ica_gc.writeToBed(OUTPUT_GC);
			} catch (IOException e) {
				String msg = String.format("Unable to write to %s", OUTPUT_GC);
				log.error(e, msg);
				throw new RuntimeException(msg, e);
			}
		}
	}
	public static void main(String[] argv) {
        System.exit(new ComputeCoverage().instanceMain(argv));
    }
}
