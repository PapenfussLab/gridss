package gridss;

import au.edu.wehi.idsv.ReadPairConcordanceMethod;
import gridss.analysis.CollectGridssMetrics;
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
    @Argument(shortName="MO", doc="Output file containing SV metrics", optional=true)
    public File REGION_BED;
    @Argument(doc = "Extract reads whose mate maps to an export region. " +
            "If the MC tag is not present, only the starting alignment position of the mate is considered. " +
            "When determining whether the mate maps to an export region " +
            "only the primary alignment of that mate is considered. Secondary " +
            "and supplementary alignments are ignored.")
    public boolean EXTRACT_MATES = true;
    @Argument(doc = "Extract all records for reads that have a chimeric alignment mapping to an export region")
    public boolean EXTRACT_SPLITS = true;
    @Argument(doc = "Number of bases surrounding each export region to include in the index query. " +
            "Setting this to a value slightly greater than the 99.99% fragment size distribution will reduce the number of random access" +
            "IO requests made. " +
            "This parameter is not used if a linear scan is performed.")
    public int REGION_PADDING_SIZE = 2000;
    @Argument(doc = "File to write the output reads to.")
    public File READ_OUTPUT;
    @Argument(doc="Number of worker threads to spawn. Defaults to number of cores available."
            + " Note that I/O threads are not included in this worker thread count so CPU usage can be higher than the number of worker thread.",
            shortName="THREADS")
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
    @Override
    public void setProgramsToRun(Collection<ProgramInterface> programsToRun) {
    	// Inject SV read extraction
    	programsToRun.add(createExtractFullReads());
    	super.setProgramsToRun(programsToRun);
    }
}
