/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecord.SAMTagAndValue;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.analysis.SinglePassSamProgram;

import java.io.File;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;

@CommandLineProgramProperties(
		summary = "Reads a SAM or BAM file and writes a file containing metrics about " +
                "the presence of SAM tags",
        oneLineSummary = "Writes SAM tag metrics for a SAM or BAM file",
        programGroup = gridss.cmdline.programgroups.Metrics.class
)
public class CollectTagMetrics extends SinglePassSamProgram {
	public static final String METRICS_SUFFIX = ".tag_metrics";
	
	@Argument(doc="If true, also include reads marked as duplicates.")
	public boolean INCLUDE_DUPLICATES = false;
	
	private Map<String, TagSummaryMetrics> tags = new HashMap<>();

    /** Required main method. */
    public static void main(final String[] args) {
        System.exit(new CollectTagMetrics().instanceMain(args));
    }

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        IOUtil.assertFileIsWritable(OUTPUT);
        tags = new HashMap<>();
    }

    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
    	if (rec.getDuplicateReadFlag() && !INCLUDE_DUPLICATES) return;
    	for (SAMTagAndValue attr : rec.getAttributes()) {
    		String tag = attr.tag;
    		TagSummaryMetrics metric = tags.get(tag);
    		if (metric == null) {
    			metric = new TagSummaryMetrics();
    			metric.TAG = tag;
    			metric.COUNT = 0;
    			tags.put(tag, metric);
    		}
    		metric.COUNT++;
    		
    	}
    }
    
    @Override
    protected void finish() {
        final MetricsFile<TagSummaryMetrics, Integer> metrics = getMetricsFile();
        tags.values().stream()
        	.sorted(Comparator.comparing(m -> m.TAG))
        	.forEach(metric -> { metrics.addMetric(metric); });
        metrics.write(OUTPUT);
    }
}
