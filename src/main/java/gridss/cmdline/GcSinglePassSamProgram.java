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
	@Argument(doc="Minimum concordant read pair fragment size if using the fixed method of calculation", optional=true)
	public Integer READ_PAIR_CONCORDANCE_MIN_FRAGMENT_SIZE = null;
	@Argument(doc="Maximum concordant read pair fragment size if using the fixed method of calculation", optional=true)
	public Integer READ_PAIR_CONCORDANCE_MAX_FRAGMENT_SIZE = null;
	@Argument(doc = "Percent (0.0-1.0) of read pairs considered concordant if using the library distribution to determine concordance.", optional=true)
	public Double READ_PAIR_CONCORDANT_PERCENT = null;
	@Argument(doc="Picard tools insert size distribution metrics txt file. Required if using the library distribution to determine concordance.", optional=true)
	public File INSERT_SIZE_METRICS = null;
	public String[] customCommandLineValidation_ProcessStructuralVariantReadsCommandLineProgram() {
		if ((READ_PAIR_CONCORDANCE_MAX_FRAGMENT_SIZE != null && READ_PAIR_CONCORDANCE_MIN_FRAGMENT_SIZE == null) ||
				READ_PAIR_CONCORDANCE_MAX_FRAGMENT_SIZE == null && READ_PAIR_CONCORDANCE_MIN_FRAGMENT_SIZE != null) {
			return new String[] { "READ_PAIR_CONCORDANCE_MIN_FRAGMENT_SIZE and READ_PAIR_CONCORDANCE_MAX_FRAGMENT_SIZE must both be specified." };
		} else if (READ_PAIR_CONCORDANT_PERCENT != null && INSERT_SIZE_METRICS == null) {
			return new String[] { "INSERT_SIZE_METRICS must be specified if READ_PAIR_CONCORDANT_PERCENT is specified." };
		}
		return super.customCommandLineValidation();
	}
	private ReadPairConcordanceCalculator rpcc = null;
	public ReadPairConcordanceCalculator getReadPairConcordanceCalculator() {
		if (rpcc == null) {
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
					READ_PAIR_CONCORDANCE_MIN_FRAGMENT_SIZE == null ? 0 : READ_PAIR_CONCORDANCE_MIN_FRAGMENT_SIZE,
					READ_PAIR_CONCORDANCE_MAX_FRAGMENT_SIZE == null ? 0 : READ_PAIR_CONCORDANCE_MAX_FRAGMENT_SIZE,
					READ_PAIR_CONCORDANT_PERCENT,
					isd,
					null);
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
