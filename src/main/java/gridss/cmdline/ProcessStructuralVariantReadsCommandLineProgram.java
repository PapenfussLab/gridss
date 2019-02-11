package gridss.cmdline;

import java.io.File;

import org.broadinstitute.barclay.argparser.Argument;

import au.edu.wehi.idsv.ReadPairConcordanceCalculator;
import au.edu.wehi.idsv.ReadPairConcordanceMethod;
import gridss.analysis.InsertSizeDistribution;
import htsjdk.samtools.util.Log;

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
    @Argument(doc="Include unmapped reads", optional=true)
    public boolean UNMAPPED_READS = true;
    @Argument(doc="If true, also include reads marked as duplicates.")
	public boolean INCLUDE_DUPLICATES = false;
	//#endregion SV read args

    @Override
	public String[] customCommandLineValidation() {
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
    public void copyInput(ProcessStructuralVariantReadsCommandLineProgram to) {
    	super.copyInput(to);
    	to.MIN_INDEL_SIZE = MIN_INDEL_SIZE;
    	to.MIN_CLIP_LENGTH = MIN_CLIP_LENGTH;
    	to.CLIPPED = CLIPPED;
    	to.INDELS = INDELS;
    	to.SPLIT = SPLIT;
    	to.SINGLE_MAPPED_PAIRED = SINGLE_MAPPED_PAIRED;
    	to.DISCORDANT_READ_PAIRS = DISCORDANT_READ_PAIRS;
    	to.READ_PAIR_CONCORDANCE_METHOD = READ_PAIR_CONCORDANCE_METHOD;
    	to.FIXED_READ_PAIR_CONCORDANCE_MIN_FRAGMENT_SIZE = FIXED_READ_PAIR_CONCORDANCE_MIN_FRAGMENT_SIZE;
    	to.FIXED_READ_PAIR_CONCORDANCE_MAX_FRAGMENT_SIZE = FIXED_READ_PAIR_CONCORDANCE_MAX_FRAGMENT_SIZE;
    	to.READ_PAIR_CONCORDANT_PERCENT = READ_PAIR_CONCORDANT_PERCENT;
    	to.INSERT_SIZE_METRICS = INSERT_SIZE_METRICS;
    	to.UNMAPPED_READS = UNMAPPED_READS;
    }
    @Override
	public boolean referenceRequired() { return false; }
}
