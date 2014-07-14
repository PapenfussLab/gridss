package au.edu.wehi.idsv;

import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.IOException;

import picard.analysis.InsertSizeMetrics;
import picard.cmdline.Usage;

import com.google.common.collect.Iterables;

/**
 * Extracts structural variation evidence and assembles breakends
 * @author Daniel Cameron
 *
 */
public class ExtractAndAssemble extends CommandLineProgram {
	private static final String PROGRAM_VERSION = "0.1";

    // The following attributes define the command-line arguments
    @Usage
    public String USAGE = getStandardUsagePreamble() + "Extracts structural variation evidence and assembles breakend FASTQs requiring realignment back to reference" + PROGRAM_VERSION;
    private Log log = Log.getInstance(ExtractAndAssemble.class);
    @Override
	protected int doWork() {
    	try {
	    	ensureDictionariesMatch();
	    	log.info("Extracting evidence");
	    	for (File f : Iterables.concat(INPUT, INPUT_TUMOUR)) {
	    		SAMEvidenceSource sref = new SAMEvidenceSource(getContext(), f);
	    		sref.ensureEvidenceExtracted();
	    	}
	    	log.info("Evidence extraction complete.");
	    	log.info("Starting assembly");
	    	throw new RuntimeException("NYI: GenerateDirectedBreakpoints");
    	} catch (IOException e) {
    		log.error(e);
    		throw new RuntimeException(e);
    	}
    }
	public static void main(String[] argv) {
        System.exit(new ExtractAndAssemble().instanceMain(argv));
    }
}
