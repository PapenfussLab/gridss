package au.edu.wehi.socrates;

import java.io.File;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.filter.AggregateFilter;
import net.sf.picard.filter.DuplicateReadFilter;
import net.sf.picard.filter.FailsVendorReadQualityFilter;
import net.sf.picard.filter.FilteringIterator;
import net.sf.picard.filter.SamRecordFilter;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.picard.util.ProgressLogger;
import net.sf.samtools.SAMException;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.util.CloseableIterator;

public class ExtractEvidence extends CommandLineProgram {

    private static final String PROGRAM_VERSION = "0.1";

    // The following attributes define the command-line arguments
    @Usage
    public String USAGE = getStandardUsagePreamble() + "Extract reads supporting structural variation." + PROGRAM_VERSION;

    @Option(doc="The BAM file to process.",
    		shortName=StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;
    @Option(doc = "Coordinate sorted output file for soft clipped reads",
            optional = false,
            shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File SC_OUTPUT;
    @Option(doc = "Coordinate sorted output file for discordantly paired reads",
            optional = false,
            shortName = "DPO")
    public File DP_OUTPUT = null;
    @Option(doc = "Coordinate sorted output file for reads pairs with only a single read mapped (open ended anchor)",
            optional = false,
            shortName = "OEAO")
    public File OEA_OUTPUT = null;
    @Option(doc = "Extract supporting reads into a separate file for each chromosome.",
            optional = true,
            shortName = "BYCHR")
    public boolean PER_CHR = true;
    @Option(doc = "Extract supporting reads into a separate file for each read group.",
            optional = true,
            shortName = "BYRG")
    public boolean PER_READ_GROUP = false;
    @Option(doc = "DP and OEA reads sorted by coordinate of mapped mate read.",
            optional = true,
            shortName = "MCO")
    public File MATE_COORDINATE_OUTPUT = null;
    
    private Log log = Log.getInstance(ExtractEvidence.class);
    @Override
	protected int doWork() {
    	// Sort output
        return 0;
    }
	public static void main(String[] argv) {
        System.exit(new ExtractEvidence().instanceMain(argv));
    }
}
