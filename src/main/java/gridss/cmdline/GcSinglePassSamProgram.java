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

package gridss.cmdline;

import au.edu.wehi.idsv.ReadPairConcordanceCalculator;
import au.edu.wehi.idsv.ReadPairConcordanceMethod;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.picard.TwoBitBufferedReferenceSequenceFile;
import gridss.analysis.InsertSizeDistribution;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import picard.analysis.SinglePassSamProgram;

import java.io.File;
import java.io.FileNotFoundException;

public abstract class GcSinglePassSamProgram extends SinglePassSamProgram {
	private static final Log log = Log.getInstance(GcSinglePassSamProgram.class);
    @Argument(doc="Fragment size to use when inferring fragment GC content of single-end or discordantly paired reads", optional=false)
    public int UNPAIRED_FRAGMENT_SIZE;
    // --------- start chunk from ProcessStructuralVariantReadsCommandLineProgram ---------
    @Argument(doc="Method of calculating read pair concordance. Valid values are SAM_FLAG, PERCENTAGE, and FIXED", optional=true)
    public ReadPairConcordanceMethod READ_PAIR_CONCORDANCE_METHOD = ReadPairConcordanceMethod.SAM_FLAG;
    @Argument(doc="Minimum concordant read pair fragment size if using the FIXED method of calculation", optional=true)
    public int FIXED_READ_PAIR_CONCORDANCE_MIN_FRAGMENT_SIZE = 0;
    @Argument(doc="Maximum concordant read pair fragment size if using the FIXED method of calculation", optional=true)
    public int FIXED_READ_PAIR_CONCORDANCE_MAX_FRAGMENT_SIZE = 0;
    @Argument(doc = "Percent (0.0-1.0) of read pairs considered concorant if using the PERCENTAGE method of calculation.", optional=true)
    public Float READ_PAIR_CONCORDANT_PERCENT = 0.995f;
    @Argument(doc="Picard tools insert size distribution metrics txt file. Required if using the PERCENTAGE read pair concordance calculation method.", optional=true)
    public File INSERT_SIZE_METRICS = null;
	protected String[] customCommandLineValidation_ProcessStructuralVariantReadsCommandLineProgram() {
		if (READ_PAIR_CONCORDANCE_METHOD == ReadPairConcordanceMethod.PERCENTAGE && INSERT_SIZE_METRICS == null) {
			return new String[] { "INSERT_SIZE_METRICS is required when using percentage based read pair concordance" };
		}
		return super.customCommandLineValidation();
	}
    private ReadPairConcordanceCalculator rpcc = null;
    private boolean rpccInitialised = false;
    public ReadPairConcordanceCalculator getReadPairConcordanceCalculator() {
    	if (!rpccInitialised) {
    		// Read metrics file
    		InsertSizeDistribution isd = null;
    		if (INSERT_SIZE_METRICS != null) {
    			if (!INSERT_SIZE_METRICS.exists()) {
    				log.warn("Missing " + INSERT_SIZE_METRICS);
    			} else {
    				isd = InsertSizeDistribution.create(INSERT_SIZE_METRICS);
    			}
    		}
    		rpcc = ReadPairConcordanceCalculator.create(
    				READ_PAIR_CONCORDANCE_METHOD,
    				FIXED_READ_PAIR_CONCORDANCE_MIN_FRAGMENT_SIZE,
    				FIXED_READ_PAIR_CONCORDANCE_MAX_FRAGMENT_SIZE,
    				READ_PAIR_CONCORDANT_PERCENT,
    				isd,
    				null);
    		rpccInitialised = true;
    	}
    	return rpcc;
    }
    // --------- end chunk from ProcessStructuralVariantReadsCommandLineProgram ---------
    // --------- start chunk from ReferenceCommandLineProgram ---------
    @Argument(doc="If true, also include reads marked as duplicates.")
    public boolean INCLUDE_DUPLICATES = false;
	private ReferenceLookup reference;
	public ReferenceLookup getReference() {
		IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
		if (reference == null) {
			try {
				reference = new TwoBitBufferedReferenceSequenceFile(new IndexedFastaSequenceFile(REFERENCE_SEQUENCE));
			} catch (FileNotFoundException e) {
				String msg = String.format("Missing reference genome %s", REFERENCE_SEQUENCE);
				log.error(msg);
				throw new RuntimeException(msg);
			}
		}
		return reference;
	}
	public void setReference(ReferenceLookup ref) {
		this.reference = ref;
	}
	@Override
	protected String[] customCommandLineValidation() {
		String[] val = referenceCustomCommandLineValidation();
		if (val != null) {
			return val;
		}
		val = customCommandLineValidation_ProcessStructuralVariantReadsCommandLineProgram();
		if (val != null) {
			return val;
		}
		return super.customCommandLineValidation();
	}
	public String[] referenceCustomCommandLineValidation() {
		if (referenceRequired()) {
			if (REFERENCE_SEQUENCE == null) {
	            return new String[]{"Must have a non-null REFERENCE_SEQUENCE"};
	        }
		}
		return null;
	}
	public boolean referenceRequired() { return true; }
    // --------- end chunk from ReferenceCommandLineProgram ---------
}
