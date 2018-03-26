package gridss.analysis;

import java.io.File;
import java.util.Set;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import picard.PicardException;
import picard.analysis.MetricAccumulationLevel;
import picard.analysis.SinglePassSamProgram;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;
import picard.util.RExecutor;

/**
 * Command line program to read non-duplicate insert sizes, create a Histogram
 * and report distribution statistics.
 *
 * @author Daniel Cameron
 */
@CommandLineProgramProperties(
		summary = "Reads a SAM or BAM file and writes a file containing metrics about " +
                "the statistical distribution of read mapping qualities (excluding duplicates) " +
                "and generates a Histogram plot.",
        oneLineSummary = "Writes mapq distribution metrics for a SAM or BAM file",
        programGroup = DiagnosticsAndQCProgramGroup.class
)
public class CollectMapqMetrics extends SinglePassSamProgram {
	public static final String METRICS_SUFFIX = ".mapq_metrics";
	public static final String HISTOGRAM_SUFFIX = ".mapq_histogram.pdf";
    private static final String Histogram_R_SCRIPT = "gridss/analysis/mapqHistogram.R";
    
    @Argument(doc="If true, also include reads marked as duplicates.")
    public boolean INCLUDE_DUPLICATES = false;

    @Argument(shortName="H", doc="File to write insert size Histogram chart to.")
    public File Histogram_FILE = null;

    @Argument(shortName="LEVEL", doc="The level(s) at which to accumulate metrics.  ")
    private Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS);

    // Calculates MapqMetrics for all METRIC_ACCUMULATION_LEVELs provided
    private MapqMetricsCollector multiCollector;

    /** Required main method implementation. */
    public static void main(final String[] argv) {
        new CollectMapqMetrics().instanceMainWithExit(argv);
    }

    /**
     * Put any custom command-line validation in an override of this method.
     * clp is initialized at this point and can be used to print usage and access argv.
     * Any options set by command-line parser can be validated.
     *
     * @return null if command line is valid.  If command line is invalid, returns an array of error message
     *         to be written to the appropriate place.
     */
    @Override
    protected String[] customCommandLineValidation() {
         return super.customCommandLineValidation();
    }

    @Override protected void setup(final SAMFileHeader header, final File samFile) {
        IOUtil.assertFileIsWritable(OUTPUT);
        if (Histogram_FILE != null) {
        	IOUtil.assertFileIsWritable(Histogram_FILE);
        }

        //Delegate actual collection to MapqMetricCollector
        multiCollector = new MapqMetricsCollector(METRIC_ACCUMULATION_LEVEL, header.getReadGroups());
    }

    @Override protected void acceptRead(final SAMRecord record, final ReferenceSequence ref) {
    	if (record.getDuplicateReadFlag() && !INCLUDE_DUPLICATES) return;
        multiCollector.acceptRecord(record, ref);
    }

    @Override protected void finish() {
        multiCollector.finish();

        final MetricsFile<MapqMetrics, Integer> file = getMetricsFile();
        multiCollector.addAllLevelsToFile(file);

        file.write(OUTPUT);

        if (Histogram_FILE != null) {
        	final int rResult = RExecutor.executeFromClasspath(
                Histogram_R_SCRIPT,
                OUTPUT.getAbsolutePath(),
                Histogram_FILE.getAbsolutePath(),
                INPUT.getName());
	        if (rResult != 0) {
	            throw new PicardException("R script " + Histogram_R_SCRIPT + " failed with return code " + rResult);
	        }
        }        
    }
}
