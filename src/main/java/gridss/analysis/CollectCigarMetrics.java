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

import au.edu.wehi.idsv.sam.CigarUtil;
import com.google.common.collect.Iterables;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.analysis.SinglePassSamProgram;

import java.io.File;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.List;

@CommandLineProgramProperties(
        summary = "Reads a SAM or BAM file and writes a file containing metrics about " +
                "the statistical distribution of alignment CIGARs.",
        oneLineSummary = "Writes CIGAR distribution metrics for a SAM or BAM file",
        programGroup = gridss.cmdline.programgroups.Metrics.class
)
public class CollectCigarMetrics extends SinglePassSamProgram {
	public static final String METRICS_SUFFIX = ".cigar_metrics";
	
	@Argument(shortName="Z", doc="If set to true include a zero length operator for each operator not included in the alignment CIGAR.")
    public boolean INCLUDE_OMITTED_OPERATORS = true;
	
	@Argument(doc="If true, also include reads marked as duplicates.")
    public boolean INCLUDE_DUPLICATES = false;
	
    private EnumMap<CigarOperator, List<CigarDetailMetrics>> cigar;    

    /** Required main method. */
    public static void main(final String[] args) {
        System.exit(new CollectCigarMetrics().instanceMain(args));
    }

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        IOUtil.assertFileIsWritable(OUTPUT);
        HashMap<CigarOperator, List<CigarDetailMetrics>> hm = new HashMap<CigarOperator, List<CigarDetailMetrics>>();
        for (CigarOperator op : CigarOperator.values()) {
			hm.put(op, new ArrayList<CigarDetailMetrics>());
		}
        cigar = new EnumMap<>(hm);
    }

    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
    	// Skip unwanted records
    	if (rec.getReadUnmappedFlag()) return;
    	if (rec.getCigar() == null) return;
    	if (rec.getDuplicateReadFlag() && !INCLUDE_DUPLICATES) return;
    	List<CigarElement> list = rec.getCigar().getCigarElements();
    	if (list == null || list.size() == 0) return;
    	for (CigarElement ce : list) {
    		acceptCigarElement(ce);
    	}
    	if (INCLUDE_OMITTED_OPERATORS) {
	    	for (CigarOperator op : CigarOperator.values()) {
	    		switch (op) {
	    			case S:
	    				if (CigarUtil.getStartSoftClipLength(list) == 0) {
	    					acceptCigarElement(new CigarElement(0, CigarOperator.S));
	    				}
	    				if (CigarUtil.getEndSoftClipLength(list) == 0) {
	    					acceptCigarElement(new CigarElement(0, CigarOperator.S));
	    				}
	    				break;
	    			case H:
	    				if (list.get(0).getOperator() != CigarOperator.H) {
	    					acceptCigarElement(new CigarElement(0, CigarOperator.H));
	    				}
	    				if (list.get(list.size() - 1).getOperator() != CigarOperator.H) {
	    					acceptCigarElement(new CigarElement(0, CigarOperator.H));
	    				}
	    				break;
	    			default:
	    				if (!Iterables.any(list, ce -> ce.getOperator() == op)) {
	    					acceptCigarElement(new CigarElement(0, op));
	    				}
	    				break;
	    		}
	    	}
    	}
    }
    
    private void acceptCigarElement(CigarElement ce) {
    	List<CigarDetailMetrics> list = cigar.get(ce.getOperator());
    	int length = ce.getLength();
    	while (list.size() <= length) {
    		CigarDetailMetrics cdm = new CigarDetailMetrics();
    		cdm.LENGTH = list.size();
    		cdm.OPERATOR = (char)CigarOperator.enumToCharacter(ce.getOperator());
    		cdm.COUNT = 0;
    		list.add(cdm);
    	}
    	list.get(ce.getLength()).COUNT++;
	}
    
    @Override
    protected void finish() {
    	// TODO: build histograms?
        final MetricsFile<CigarDetailMetrics, Integer> metrics = getMetricsFile();
        cigar.values().stream().flatMap(c -> c.stream()).forEach(metric -> {
        	metrics.addMetric(metric);
		});
        metrics.write(OUTPUT);
    }
}
