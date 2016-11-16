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

package gridss.analysis.directed;

import java.util.List;
import java.util.Set;

import gridss.analysis.MapqMetrics;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Histogram;
import picard.analysis.MetricAccumulationLevel;
import picard.metrics.MultiLevelCollector;
import picard.metrics.PerUnitMetricCollector;

/**
 * Collects InserSizeMetrics on the specified accumulationLevels using
 */
public class MapqMetricsCollector extends MultiLevelCollector<MapqMetrics, Integer, Integer> {

    public MapqMetricsCollector(final Set<MetricAccumulationLevel> accumulationLevels, final List<SAMReadGroupRecord> samRgRecords) {
        setup(accumulationLevels, samRgRecords);
    }

    @Override
    protected Integer makeArg(SAMRecord samRecord, ReferenceSequence refSeq) {
    	return samRecord.getMappingQuality();
    }

    /** Make an InsertSizeCollector with the given arguments */
    @Override
    protected PerUnitMetricCollector<MapqMetrics, Integer, Integer> makeChildCollector(final String sample, final String library, final String readGroup) {
        return new PerUnitMapqMetricsCollector(sample, library, readGroup);
    }

    @Override
    public void acceptRecord(final SAMRecord record, final ReferenceSequence refSeq) {
        if (record.getReadUnmappedFlag()) {
            return;
        }
        super.acceptRecord(record, refSeq);
    }

    /** A Collector for individual InsertSizeMetrics for a given SAMPLE or SAMPLE/LIBRARY or SAMPLE/LIBRARY/READ_GROUP (depending on aggregation levels) */
    public class PerUnitMapqMetricsCollector implements PerUnitMetricCollector<MapqMetrics, Integer, Integer> {
        final Histogram<Integer> histogram;
        final String sample;
        final String library;
        final String readGroup;

        public PerUnitMapqMetricsCollector(final String sample, final String library, final String readGroup) {
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
            histogram = new Histogram<Integer>("mapq", prefix);
        }

        public void acceptRecord(final Integer args) {
        	histogram.increment(args);
        }

        public void finish() { }

        public void addMetricsToFile(final MetricsFile<MapqMetrics,Integer> file) {
        	final MapqMetrics metrics = new MapqMetrics();
            metrics.SAMPLE             = this.sample;
            metrics.LIBRARY            = this.library;
            metrics.READ_GROUP         = this.readGroup;
            metrics.MAPPED_READS       = (long) histogram.getCount();
            metrics.MIN_MAPQ           = SAMRecord.UNKNOWN_MAPPING_QUALITY;
            metrics.MAX_MAPQ           = 0;
            for (int key : histogram.keySet()) {
            	if (key != SAMRecord.UNKNOWN_MAPPING_QUALITY) {
            		if (histogram.get(key).getValue() > 0) {
            			metrics.MIN_MAPQ = Math.min(metrics.MIN_MAPQ, key);
            			metrics.MAX_MAPQ = Math.max(metrics.MAX_MAPQ, key);
            		}
            	}
            }
            if (histogram.get(SAMRecord.NO_MAPPING_QUALITY) != null) {
            	metrics.ZERO_MAPQ      = (long) histogram.get(SAMRecord.NO_MAPPING_QUALITY).getValue();
            }
            if (histogram.get(SAMRecord.UNKNOWN_MAPPING_QUALITY) != null) {
                metrics.UNKNOWN_MAPQ   = (long) histogram.get(SAMRecord.UNKNOWN_MAPPING_QUALITY).getValue();
            }
            file.addHistogram(histogram);
            file.addMetric(metrics);       
        }
    }
}

