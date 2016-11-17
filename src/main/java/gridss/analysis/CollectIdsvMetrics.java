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

import java.io.File;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamPairUtil.PairOrientation;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.IOUtil;
import picard.analysis.SinglePassSamProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.programgroups.Metrics;

@CommandLineProgramProperties(
        usage = "Reads a SAM or BAM file and writes a file containing metrics " +
                "used by idsv.",
        usageShort = "Writes idsv metrics for a SAM or BAM file",
        programGroup = Metrics.class
)
public class CollectIdsvMetrics extends SinglePassSamProgram {
	
    private IdsvMetrics idsv;    

    /** Required main method. */
    public static void main(final String[] args) {
        System.exit(new CollectIdsvMetrics().instanceMain(args));
    }

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        IOUtil.assertFileIsWritable(OUTPUT);
        idsv = new IdsvMetrics();
    }

    @Override
    protected void acceptRead(final SAMRecord record, final ReferenceSequence ref) {
    	idsv.MAX_READ_LENGTH = Math.max(idsv.MAX_READ_LENGTH, record.getReadLength());
    	if (!record.getReadUnmappedFlag()) {
    		idsv.MAX_READ_MAPPED_LENGTH = Math.max(idsv.MAX_READ_MAPPED_LENGTH, record.getAlignmentEnd() - record.getAlignmentStart() + 1);
    	}
    	if (record.getReadPairedFlag()) {
    		if (record.getProperPairFlag()) {
	    		int fragmentSize = SAMRecordUtil.estimateFragmentSize(record, PairOrientation.FR);
	    		fragmentSize = Math.abs(fragmentSize);
	    		if (idsv.MAX_PROPER_PAIR_FRAGMENT_LENGTH == null) {
	    			idsv.MAX_PROPER_PAIR_FRAGMENT_LENGTH = fragmentSize;
	    		} else {
	    			idsv.MAX_PROPER_PAIR_FRAGMENT_LENGTH = Math.max(idsv.MAX_PROPER_PAIR_FRAGMENT_LENGTH, Math.abs(fragmentSize));
	    		}
	    		if (idsv.MIN_PROPER_PAIR_FRAGMENT_LENGTH == null) {
	    			idsv.MIN_PROPER_PAIR_FRAGMENT_LENGTH = fragmentSize;
	    		} else {
	    			idsv.MIN_PROPER_PAIR_FRAGMENT_LENGTH = Math.min(idsv.MIN_PROPER_PAIR_FRAGMENT_LENGTH, Math.abs(fragmentSize));
	    		}
    		}
    		if (record.getFirstOfPairFlag()) {
    			idsv.READ_PAIRS++;
    			if (record.getReadUnmappedFlag() && record.getMateUnmappedFlag()) {
    				idsv.READ_PAIRS_ZERO_MAPPED++;
    			} else if (!record.getReadUnmappedFlag() && !record.getMateUnmappedFlag()) {
    				idsv.READ_PAIRS_BOTH_MAPPED++;
    			} else {
    				idsv.READ_PAIRS_ONE_MAPPED++;
    			}
    		}
    	}
    	idsv.READS++;
    	if (!record.getReadUnmappedFlag()) {
    		idsv.MAPPED_READS++;
    	}
    }
    
    @Override
    protected void finish() {
        final MetricsFile<IdsvMetrics, Integer> metricsFile = getMetricsFile();
        metricsFile.addMetric(idsv);
        metricsFile.write(OUTPUT);
    }
}
