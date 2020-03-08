package gridss.cmdline;

import au.edu.wehi.idsv.ReadPairConcordanceCalculator;
import gridss.analysis.InsertSizeDistribution;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;

import java.io.File;

/**
 * Base class used to transform a VCF breakpoint call set given the full evidence available.
 * 
 * 
 * @author Daniel Cameron
 *
 */
public abstract class ProcessStructuralVariantReadsCommandLineProgram extends ByReadNameSinglePassSamProgram {
	private static final Log log = Log.getInstance(ProcessStructuralVariantReadsCommandLineProgram.class);
	//#region SV read args
	@Argument(doc="Minimum indel size", optional=true)
    public int MIN_INDEL_SIZE = 1;
    @Argument(doc="Minimum bases clipped", optional=true)
    public int MIN_CLIP_LENGTH = 1;
    @Argument(doc="Include hard and soft clipped reads in output", optional=true)
    public boolean CLIPPED = true;
    @Argument(doc="Include reads containing indels in output", optional=true)
    public boolean INDELS = true;
    @Argument(doc="Include split reads in output", optional=true)
    public boolean SPLIT = true;
    @Argument(doc="Include read pairs in which only one of the read is aligned to the reference.", optional=true)
    public boolean SINGLE_MAPPED_PAIRED = true;
    @Argument(doc="Include read pairs that align do not align in the expected orientation within the expected fragment size distribution.", optional=true)
    public boolean DISCORDANT_READ_PAIRS = true;
    @Argument(doc="Minimum concordant read pair fragment size if using the fixed method of calculation", optional=true)
    public Integer READ_PAIR_CONCORDANCE_MIN_FRAGMENT_SIZE = null;
    @Argument(doc="Maximum concordant read pair fragment size if using the fixed method of calculation", optional=true)
    public Integer READ_PAIR_CONCORDANCE_MAX_FRAGMENT_SIZE = null;
    @Argument(doc = "Percent (0.0-1.0) of read pairs considered concordant if using the library distribution to determine concordance.", optional=true)
    public Double READ_PAIR_CONCORDANT_PERCENT = null;
    @Argument(doc="Picard tools insert size distribution metrics txt file. Required if using the library distribution to determine concordance.", optional=true)
    public File INSERT_SIZE_METRICS = null;
    @Argument(doc="Include unmapped reads", optional=true)
    public boolean UNMAPPED_READS = true;
    @Argument(doc="If true, also include reads marked as duplicates.")
	public boolean INCLUDE_DUPLICATES = false;
	//#endregion SV read args

    @Override
	public String[] customCommandLineValidation() {
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
    public void copyInput(ProcessStructuralVariantReadsCommandLineProgram to) {
    	super.copyInput(to);
    	to.MIN_INDEL_SIZE = MIN_INDEL_SIZE;
    	to.MIN_CLIP_LENGTH = MIN_CLIP_LENGTH;
    	to.CLIPPED = CLIPPED;
    	to.INDELS = INDELS;
    	to.SPLIT = SPLIT;
    	to.SINGLE_MAPPED_PAIRED = SINGLE_MAPPED_PAIRED;
    	to.DISCORDANT_READ_PAIRS = DISCORDANT_READ_PAIRS;
    	to.READ_PAIR_CONCORDANCE_MIN_FRAGMENT_SIZE = READ_PAIR_CONCORDANCE_MIN_FRAGMENT_SIZE;
    	to.READ_PAIR_CONCORDANCE_MAX_FRAGMENT_SIZE = READ_PAIR_CONCORDANCE_MAX_FRAGMENT_SIZE;
    	to.READ_PAIR_CONCORDANT_PERCENT = READ_PAIR_CONCORDANT_PERCENT;
    	to.INSERT_SIZE_METRICS = INSERT_SIZE_METRICS;
    	to.UNMAPPED_READS = UNMAPPED_READS;
    	to.rpcc = rpcc;
    }
    @Override
	public boolean referenceRequired() { return false; }
}
