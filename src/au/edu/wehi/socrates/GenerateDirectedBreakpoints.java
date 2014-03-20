package au.edu.wehi.socrates;

import java.io.File;

import org.apache.commons.cli.OptionBuilder;

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

public class GenerateDirectedBreakpoints extends CommandLineProgram {

    private static final String PROGRAM_VERSION = "0.1";

    // The following attributes define the command-line arguments
    @Usage
    public String USAGE = getStandardUsagePreamble() + "Generated directed breakpoints." + PROGRAM_VERSION;

    @Option(doc = "Coordinate sorted input file containing soft clipped reads",
            optional = false,
            shortName = "SCI")
    public File SC_INPUT;
    @Option(doc = "Coordinate sorted input file containing discordantly paired reads",
            optional = false,
            shortName = "DPI")
    public File DP_INPUT = null;
    @Option(doc = "Coordinate sorted input file containing reads pairs with only a single read mapped (open ended anchor)",
            optional = false,
            shortName = "OEAI")
    public File OEA_INTPUT = null;
    @Option(doc = "DP and OEA reads sorted by coordinate of mapped mate read.",
            optional = false,
            shortName = "MCI")
    public File MATE_COORDINATE_INPUT = null;    
    @Option(doc = "Directed single-ended breakpoints. A placeholder contig is output as the breakpoint partner.",
            optional = false,
            shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File VCF_OUTPUT;
    @Option(doc = "FASTQ of reference strand sequence of breakpoint excluding anchor. These sequences are used to align breakpoints",
            optional = false,
            shortName = "FQ")
    public File FASTQ_OUTPUT = null;

    @Option(doc = "Minimum alignment mapq", optional=true)
    public int MIN_MAPQ = 5;
    @Option(doc = "Length threshold of long soft-clip", optional=true)
    public int LONG_SC_LEN = 25;
    @Option(doc = "Minimum alignment percent identity to reference. Takes values in the range 0-100.", optional=true)
    public int MIN_PERCENT_IDENTITY = 95;
    @Option(doc = "Minimum average base quality score of soft clipped sequence", optional=true)
    public int MIN_LONG_SC_BASE_QUALITY = 5;
    private Log log = Log.getInstance(GenerateDirectedBreakpoints.class);
    @Override
	protected int doWork() {
    	return -1;
    	// option 1:
    	while (iter.hasNext()) {
    		process(iter.next());
    	}
    	// option 2:
    	process(iter);
    }
	public static void main(String[] argv) {
        System.exit(new GenerateDirectedBreakpoints().instanceMain(argv));
    }
}
