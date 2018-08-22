package gridss;

import java.io.File;
import java.util.Collection;
import java.util.Set;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;

import au.edu.wehi.idsv.ReadPairConcordanceMethod;
import gridss.analysis.CollectGridssMetrics;
import gridss.cmdline.CommandLineProgramHelper;
import picard.analysis.MetricAccumulationLevel;
import picard.analysis.SinglePassSamProgram;

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
public class CollectGridssMetricsAndExtractSVReads extends CollectGridssMetrics {
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
    @Argument(shortName="MO", doc="Output file containing SV metrics", optional=true)
    public File METRICS_OUTPUT;
    @Argument(shortName = "SVO", doc = "File to write the output to.")
    public File SV_OUTPUT;
    public static void main(final String[] args) {
        new CollectGridssMetricsAndExtractSVReads().instanceMainWithExit(args);
    }
    @Override
    protected String[] customCommandLineValidation() {
        return super.customCommandLineValidation();
    }
    protected ExtractSVReads getExtractSVReads() {
    	ExtractSVReads extract = new ExtractSVReads();
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
    	extract.METRICS_OUTPUT = METRICS_OUTPUT;
    	extract.OUTPUT = SV_OUTPUT;
    	extract.INPUT = this.INPUT;
    	extract.ASSUME_SORTED = true;
    	return extract;
    }
    public ProgramInterface createExtractSVReads() {
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
				final ExtractSVReads program = getExtractSVReads();
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
    	// Inject SV read extraction
    	programsToRun.add(createExtractSVReads());
    	super.setProgramsToRun(programsToRun);
    }
}
