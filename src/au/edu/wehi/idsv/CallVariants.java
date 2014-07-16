package au.edu.wehi.idsv;

import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.omg.CORBA.Environment;

import picard.cmdline.Usage;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;

/**
 * Extracts structural variation evidence and assembles breakends
 * @author Daniel Cameron
 *
 */
public class CallVariants extends CommandLineProgram {
	private static final String PROGRAM_VERSION = "0.1";

    // The following attributes define the command-line arguments
    @Usage
    public String USAGE = getStandardUsagePreamble() + "Calls variants from the extracted and assembled evidence " + PROGRAM_VERSION;
    private Log log = Log.getInstance(CallVariants.class);
    @Override
	protected int doWork() {
    	try {
	    	ensureDictionariesMatch();
	    	List<SAMEvidenceSource> samEvidence = Lists.newArrayList();
	    	for (File f : INPUT) {
	    		SAMEvidenceSource sref = new SAMEvidenceSource(getContext(), f, false);
	    		sref.ensureEvidenceExtracted();
	    		samEvidence.add(sref);
	    	}
	    	for (File f : INPUT_TUMOUR) {
	    		SAMEvidenceSource sref = new SAMEvidenceSource(getContext(), f, true);
	    		sref.ensureEvidenceExtracted();
	    		samEvidence.add(sref);
	    	}
	    	AssemblyEvidenceSource aes = new AssemblyEvidenceSource(getContext(), samEvidence, OUTPUT);
	    	List<EvidenceSource> allEvidence = Lists.newArrayList();
	    	allEvidence.add(aes);
	    	allEvidence.addAll(samEvidence);
	    	String instructions = getAlignmentInstructions(allEvidence);
	    	if (instructions != null) {
	    		log.error("Realignment is required");
	    		log.error(instructions);
		    	log.error("Realignment is required");
		    	return -1;
	    	}
	    	try (VariantCaller caller = new VariantCaller(getContext(), OUTPUT, allEvidence)) {
	    		caller.callBreakends();
	    		caller.annotateBreakpoints();
	    	}
    	} catch (IOException e) {
    		log.error(e);
    		throw new RuntimeException(e);
    	}
		return 0;
    }
	public static void main(String[] argv) {
        System.exit(new CallVariants().instanceMain(argv));
    }
}
