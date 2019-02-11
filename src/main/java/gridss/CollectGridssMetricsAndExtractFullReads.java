package gridss;

import au.edu.wehi.idsv.ReadPairConcordanceMethod;
import gridss.analysis.CollectGridssMetrics;
import gridss.analysis.CollectIdsvMetrics;
import gridss.analysis.CollectStructuralVariantReadMetrics;
import gridss.cmdline.CommandLineProgramHelper;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.analysis.MetricAccumulationLevel;
import picard.analysis.SinglePassSamProgram;

import java.io.File;
import java.util.Collection;
import java.util.Set;

/**
 * Class that is designed to instantiate and execute multiple metrics programs that extend
 * SinglePassSamProgram while making only a single pass through the SAM file and supplying
 * each program with the records as it goes.
 *
 */
@CommandLineProgramProperties(
        summary = "Merging of CollectGridssMetrics and ExtactSVReads designed for pipelines that already have insert size metrics. "
        		+ "Combining the two programs removes the unnecessary CPU overhead of parsing the input file multiple times.",
        oneLineSummary = "A \"meta-metrics\" calculating program that produces multiple metrics for the provided SAM/BAM and extracts SV reads.",
        programGroup = gridss.cmdline.programgroups.DataConversion.class
)
public class CollectGridssMetricsAndExtractFullReads extends CollectGridssMetrics {
    // #region SV command line args
    @Argument(doc="Minimum indel size", optional=true)
    public int MIN_INDEL_SIZE = 1;
    @Argument(doc="Minimum bases clipped", optional=true)
    public int MIN_CLIP_LENGTH = 1;
    @Argument(doc="Include hard and soft clipped reads in output", optional=true)
    public boolean CLIPPED = true;
    @Argument(doc="Include reads containing indels in output", optional=true)
    public boolean INDELS = true;
    @Argument(doc="Include split reads in output", optional=true)
    public boolean SPLIT = true;
    @Argument(doc="Include read pairs in which only one of the read is aligned to the reference.", optional=true)
    public boolean SINGLE_MAPPED_PAIRED = true;
    @Argument(doc="Include read pairs that align do not align in the expected orientation within the expected fragment size distribution.", optional=true)
    public boolean DISCORDANT_READ_PAIRS = true;
    @Argument(doc="Method of calculating read pair concordance. Valid values are SAM_FLAG, PERCENTAGE, and FIXED", optional=true)
    public ReadPairConcordanceMethod READ_PAIR_CONCORDANCE_METHOD = ReadPairConcordanceMethod.SAM_FLAG;
    @Argument(doc="Minimum concordant read pair fragment size if using the FIXED method of calculation", optional=true)
    public int FIXED_READ_PAIR_CONCORDANCE_MIN_FRAGMENT_SIZE = 0;
    @Argument(doc="Maximum concordant read pair fragment size if using the FIXED method of calculation", optional=true)
    public int FIXED_READ_PAIR_CONCORDANCE_MAX_FRAGMENT_SIZE = 0;
    @Argument(doc = "Percent (0.0-1.0) of read pairs considered concorant if using the PERCENTAGE method of calculation.", optional=true)
    public Float READ_PAIR_CONCORDANT_PERCENT = 0.995f;
    @Argument(doc="Picard tools insert size distribution metrics txt file. Required if using the PERCENTAGE read pair concordance calculation method.", optional=true)
    public File INSERT_SIZE_METRICS = null;
    @Argument(doc="Include unmapped reads", optional=true)
    public boolean UNMAPPED_READS = true;
    @Argument(doc="If true, also include reads marked as duplicates.")
    public boolean INCLUDE_DUPLICATES = false;
    @Argument(shortName="SVO", doc="Output file containing SV metrics", optional=true)
    public File SV_METRICS_OUTPUT;
    // #endregion
    // #region extract command line args
    @Argument(shortName="MO", doc="Output file containing SV metrics", optional=true)
    public File REGION_BED = new ExtractFullReads().REGION_BED;
    @Argument(doc = "Extract reads whose mate maps to an export region. " +
            "If the MC tag is not present, only the starting alignment position of the mate is considered. " +
            "When determining whether the mate maps to an export region " +
            "only the primary alignment of that mate is considered. Secondary " +
            "and supplementary alignments are ignored.")
    public boolean EXTRACT_MATES = new ExtractFullReads().EXTRACT_MATES;
    @Argument(doc = "Extract all records for reads that have a chimeric alignment mapping to an export region", optional=true)
    public boolean EXTRACT_SPLITS = new ExtractFullReads().EXTRACT_SPLITS;
    @Argument(doc = "Number of bases surrounding each export region to include in the index query. " +
            "Setting this to a value slightly greater than the 99.99% fragment size distribution will reduce the number of random access" +
            "IO requests made. " +
            "This parameter is not used if a linear scan is performed.", optional=true)
    public int REGION_PADDING_SIZE = new ExtractFullReads().REGION_PADDING_SIZE;
    @Argument(doc = "File to write the output reads to.", optional=true)
    public File READ_OUTPUT;
    @Argument(doc="Number of worker threads to spawn. Defaults to number of cores available."
            + " Note that I/O threads are not included in this worker thread count so CPU usage can be higher than the number of worker thread.",
            shortName="THREADS")
    // #endregion
    public int WORKER_THREADS = 1;
    public static void main(final String[] args) {
        new CollectGridssMetricsAndExtractFullReads().instanceMainWithExit(args);
    }
    @Override
    protected String[] customCommandLineValidation() {
        return super.customCommandLineValidation();
    }
    protected ExtractFullReads getExtractFullReads() {
    	ExtractFullReads extract = new ExtractFullReads();
    	CommandLineProgramHelper.copyInputs(this, extract);
    	extract.OUTPUT = READ_OUTPUT;
        extract.INPUT = this.INPUT;
        extract.REGION_BED = REGION_BED;
        extract.EXTRACT_MATES = EXTRACT_MATES;
        extract.EXTRACT_SPLITS = EXTRACT_SPLITS;
    	extract.REGION_PADDING_SIZE = REGION_PADDING_SIZE;
    	return extract;
    }
    public ProgramInterface createExtractFullReads() {
    	return new ProgramInterface() {
			@Override
			public SinglePassSamProgram makeInstance(final String outbase,
                                                     final String outext,
                                                     final File input,
                                                     final File reference,
                                                     final Set<MetricAccumulationLevel> metricAccumulationLevel,
                                                     final File dbSnp,
                                                     final File intervals,
                                                     final File refflat,
                                                     final  Set<String> ignoreSequence) {
				final ExtractFullReads program = getExtractFullReads();
                return program;
			}
			@Override
			public boolean needsReferenceSequence() {
				return false;
			}
			@Override
			public boolean supportsMetricAccumulationLevel() {
				return false;
			}
        };
    }
    protected CollectStructuralVariantReadMetrics getSVMetrics() {
        CollectStructuralVariantReadMetrics extract = new CollectStructuralVariantReadMetrics();
        CommandLineProgramHelper.copyInputs(this, extract);
        extract.MIN_INDEL_SIZE = MIN_INDEL_SIZE;
        extract.MIN_CLIP_LENGTH = MIN_CLIP_LENGTH;
        extract.CLIPPED = CLIPPED;
        extract.INDELS = INDELS;
        extract.SPLIT = SPLIT;
        extract.SINGLE_MAPPED_PAIRED = SINGLE_MAPPED_PAIRED;
        extract.DISCORDANT_READ_PAIRS = DISCORDANT_READ_PAIRS;
        extract.READ_PAIR_CONCORDANCE_METHOD = READ_PAIR_CONCORDANCE_METHOD;
        extract.FIXED_READ_PAIR_CONCORDANCE_MIN_FRAGMENT_SIZE = FIXED_READ_PAIR_CONCORDANCE_MIN_FRAGMENT_SIZE;
        extract.FIXED_READ_PAIR_CONCORDANCE_MAX_FRAGMENT_SIZE = FIXED_READ_PAIR_CONCORDANCE_MAX_FRAGMENT_SIZE;
        extract.READ_PAIR_CONCORDANT_PERCENT = READ_PAIR_CONCORDANT_PERCENT;
        extract.INSERT_SIZE_METRICS = INSERT_SIZE_METRICS;
        extract.UNMAPPED_READS = UNMAPPED_READS;
        extract.INCLUDE_DUPLICATES = INCLUDE_DUPLICATES;
        extract.OUTPUT = SV_METRICS_OUTPUT;
        extract.INPUT = this.INPUT;
        extract.ASSUME_SORTED = true;
        return extract;
    }
    public ProgramInterface createSVMetrics() {
        return new ProgramInterface() {
            @Override
            public SinglePassSamProgram makeInstance(final String outbase,
                                                     final String outext,
                                                     final File input,
                                                     final File reference,
                                                     final Set<MetricAccumulationLevel> metricAccumulationLevel,
                                                     final File dbSnp,
                                                     final File intervals,
                                                     final File refflat,
                                                     final  Set<String> ignoreSequence) {
                final CollectStructuralVariantReadMetrics program =  getSVMetrics();
                return program.asSinglePassSamProgram();
            }
            @Override
            public boolean needsReferenceSequence() {
                return false;
            }
            @Override
            public boolean supportsMetricAccumulationLevel() {
                return false;
            }
        };
    }
    @Override
    public void setProgramsToRun(Collection<ProgramInterface> programsToRun) {
    	programsToRun.add(createExtractFullReads());
        programsToRun.add(createSVMetrics());
    	super.setProgramsToRun(programsToRun);
    }
}
