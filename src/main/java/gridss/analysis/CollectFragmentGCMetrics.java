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

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;

import gridss.cmdline.GcSinglePassSamProgram;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import picard.PicardException;
import picard.analysis.MetricAccumulationLevel;
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
        programGroup = picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup.class
)
public class CollectFragmentGCMetrics extends GcSinglePassSamProgram {
	public static final String METRICS_SUFFIX = ".gc_metrics";
	public static final String HISTOGRAM_SUFFIX = ".gc_histogram.pdf";
    private static final String Histogram_R_SCRIPT = "gridss/analysis/gcHistogram.R";
    
    @Argument(shortName="H", doc="File to write insert size Histogram chart to.")
    public File Histogram_FILE = null;

    @Argument(shortName="LEVEL", doc="The level(s) at which to accumulate metrics.  ")
    private Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS);

    // Calculates Metrics for all METRIC_ACCUMULATION_LEVELs provided
    private GcMetricsCollector multiCollector;

    @Override protected void setup(final SAMFileHeader header, final File samFile) {
        IOUtil.assertFileIsWritable(OUTPUT);
        if (Histogram_FILE != null) {
        	IOUtil.assertFileIsWritable(Histogram_FILE);
        }
        //Delegate actual collection to GcMetricsCollector
        multiCollector = new GcMetricsCollector(UNPAIRED_FRAGMENT_SIZE, getReadPairConcordanceCalculator(), METRIC_ACCUMULATION_LEVEL, header.getReadGroups());
    }

    @Override protected void acceptRead(final SAMRecord record, final ReferenceSequence ref) {
    	if (record.getDuplicateReadFlag() && IGNORE_DUPLICATES) {
    		// ignore duplicates
    	} else {
    		multiCollector.acceptRecord(record, ref);
    	}
    }

    @Override protected void finish() {
        multiCollector.finish();

        final MetricsFile<GcMetrics, Integer> file = getMetricsFile();
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
    /** Required main method implementation. */
    public static void main(final String[] argv) {
        new CollectFragmentGCMetrics().instanceMainWithExit(argv);
    }
}
