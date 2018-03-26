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
import java.io.IOException;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;

import au.edu.wehi.idsv.GenomicProcessingContext;
import au.edu.wehi.idsv.LinearGenomicCoordinate;
import au.edu.wehi.idsv.PaddedLinearGenomicCoordinate;
import au.edu.wehi.idsv.SequentialCoverageThreshold;
import au.edu.wehi.idsv.bed.IntervalBed;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequence;
import picard.analysis.SinglePassSamProgram;

@CommandLineProgramProperties(
        summary = "Reads a SAM or BAM file and writes a BED file containing the regions in which coverage equals or exceeds the given threshold",
        oneLineSummary = "Reports regions at least threshold coverage.",
        programGroup = picard.cmdline.programgroups.Metrics.class
)
public class ReportThresholdCoverage extends SinglePassSamProgram {
	public static final String SUFFIX = ".coverage.blacklist.bed";
	
	@Argument(doc = "Minimum coverage to report.", optional=false)
	public int THRESHOLD_COVERAGE;
	
	private SequentialCoverageThreshold threshold;
	
    /** Required main method. */
    public static void main(final String[] args) {
        System.exit(new ReportThresholdCoverage().instanceMain(args));
    }
    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
    	SAMSequenceDictionary dictionary = header.getSequenceDictionary();
    	LinearGenomicCoordinate linear = new PaddedLinearGenomicCoordinate(dictionary, GenomicProcessingContext.LINEAR_COORDINATE_CHROMOSOME_BUFFER, true);
    	this.threshold = new SequentialCoverageThreshold(dictionary, linear, THRESHOLD_COVERAGE);
    }

    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
    	this.threshold.acceptRecord(rec);
    }
    
    @Override
    protected void finish() {
    	IntervalBed bed = this.threshold.finish();
    	try {
			bed.write(OUTPUT, INPUT.getName());
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
    }
}
