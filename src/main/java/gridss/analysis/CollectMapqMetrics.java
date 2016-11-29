/*
 * The MIT License
 *
 * Copyright (c) 2016 Daniel Cameron
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package gridss.analysis;

import java.io.File;
import java.util.Set;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import picard.PicardException;
import picard.analysis.MetricAccumulationLevel;
import picard.analysis.SinglePassSamProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Metrics;
import picard.util.RExecutor;

/**
 * Command line program to read non-duplicate insert sizes, create a Histogram
 * and report distribution statistics.
 *
 * @author Daniel Cameron
 */
@CommandLineProgramProperties(
        usage = "Reads a SAM or BAM file and writes a file containing metrics about " +
                "the statistical distribution of read mapping qualities (excluding duplicates) " +
                "and generates a Histogram plot.",
        usageShort = "Writes mapq distribution metrics for a SAM or BAM file",
        programGroup = Metrics.class
)
public class CollectMapqMetrics extends SinglePassSamProgram {
	public static final String METRICS_SUFFIX = ".mapq_metrics";
	public static final String HISTOGRAM_SUFFIX = ".mapq_histogram.pdf";
    private static final String Histogram_R_SCRIPT = "gridss/analysis/mapqHistogram.R";

    @Option(shortName="H", doc="File to write insert size Histogram chart to.")
    public File Histogram_FILE = null;

    @Option(shortName="LEVEL", doc="The level(s) at which to accumulate metrics.  ")
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
