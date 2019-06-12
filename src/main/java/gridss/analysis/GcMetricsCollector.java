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

import au.edu.wehi.idsv.ReadGcSummary;
import au.edu.wehi.idsv.ReadPairConcordanceCalculator;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Histogram;
import picard.analysis.MetricAccumulationLevel;
import picard.metrics.MultiLevelCollector;
import picard.metrics.PerUnitMetricCollector;

import java.util.List;
import java.util.Set;

/**
 * Collects InserSizeMetrics on the specified accumulationLevels using
 */
public class GcMetricsCollector extends MultiLevelCollector<GcMetrics, Integer, Integer> {
	private final int defaultFragmentSize; 
	private final ReadPairConcordanceCalculator rpcc;
    public GcMetricsCollector(final int defaultFragmentSize, final ReadPairConcordanceCalculator rpcc,
    		final Set<MetricAccumulationLevel> accumulationLevels, final List<SAMReadGroupRecord> samRgRecords) {
    	this.defaultFragmentSize = defaultFragmentSize;
    	this.rpcc = rpcc;
        setup(accumulationLevels, samRgRecords);
    }

    @Override
    protected Integer makeArg(SAMRecord samRecord, ReferenceSequence refSeq) {
    	return (int)new ReadGcSummary(samRecord, refSeq, defaultFragmentSize, rpcc).gcPercentage;
    }

    /** Make an InsertSizeCollector with the given arguments */
    @Override
    protected PerUnitMetricCollector<GcMetrics, Integer, Integer> makeChildCollector(final String sample, final String library, final String readGroup) {
        return new PerUnitGcMetricsCollector(sample, library, readGroup);
    }

    @Override
    public void acceptRecord(final SAMRecord record, final ReferenceSequence refSeq) {
        if (record.getReadUnmappedFlag()) {
            return;
        }
        super.acceptRecord(record, refSeq);
    }

    /** A Collector for individual InsertSizeMetrics for a given SAMPLE or SAMPLE/LIBRARY or SAMPLE/LIBRARY/READ_GROUP (depending on aggregation levels) */
    public class PerUnitGcMetricsCollector implements PerUnitMetricCollector<GcMetrics, Integer, Integer> {
        final Histogram<Integer> histogram;
        final String sample;
        final String library;
        final String readGroup;

        public PerUnitGcMetricsCollector(final String sample, final String library, final String readGroup) {
            this.sample = sample;
            this.library = library;
            this.readGroup = readGroup;
            String prefix = null;
            if (this.readGroup != null) {
                prefix = this.readGroup + ".";
            }
            else if (this.library != null) {
                prefix = this.library + ".";
            }
            else if (this.sample != null) {
                prefix = this.sample + ".";
            }
            else {
                prefix = "All_Reads.";
            }
            histogram = new Histogram<Integer>("gc", prefix);
        }

        public void acceptRecord(final Integer args) {
        	histogram.increment(args);
        }

        public void finish() { }

        public void addMetricsToFile(final MetricsFile<GcMetrics,Integer> file) {
        	final GcMetrics metrics = new GcMetrics();
            metrics.SAMPLE             = this.sample;
            metrics.LIBRARY            = this.library;
            metrics.READ_GROUP         = this.readGroup;
            double totalGc = 0;
            int totalReads = 0;
            for (Integer gcBin : histogram.keySet()) {
            	int count = (int)histogram.get(gcBin).getValue();
            	totalReads += count;
            	totalGc += gcBin * count;
            }
            metrics.MEAN_GC_CONTENT = totalGc;
            metrics.READ_COUNT = totalReads;
            file.addHistogram(histogram);
            file.addMetric(metrics);       
        }
    }
}

