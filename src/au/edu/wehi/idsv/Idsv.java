package au.edu.wehi.idsv;

import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import picard.cmdline.Usage;

import com.google.common.base.Predicate;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;

/**
 * Extracts structural variation evidence and assembles breakends
 * @author Daniel Cameron
 *
 */
public class Idsv extends CommandLineProgram {
	private static final String PROGRAM_VERSION = "0.1";

    // The following attributes define the command-line arguments
    @Usage
    public String USAGE = getStandardUsagePreamble() + "Calls structural variations from NGS sequencing data. " + PROGRAM_VERSION;
    private Log log = Log.getInstance(Idsv.class);
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
	    	log.debug("Evidence extraction complete.");
	    	if (!Iterables.all(samEvidence, new Predicate<SAMEvidenceSource>() {
		    		public boolean apply(SAMEvidenceSource arg) {
						return arg.isProcessingComplete();
					}})) {
	    		log.error("Unable to extract SAM evidence. Terminating early.");
	    		return -1;
	    	}

	    	AssemblyEvidenceSource aes = new AssemblyEvidenceSource(getContext(), samEvidence, OUTPUT);
	    	aes.ensureAssembled();
	    	List<EvidenceSource> allEvidence = Lists.newArrayList();
	    	allEvidence.add(aes);
	    	allEvidence.addAll(samEvidence);
	    	if (!Iterables.all(allEvidence, new Predicate<EvidenceSource>() {
		    		public boolean apply(EvidenceSource arg) {
						return arg.isProcessingComplete();
					}})) {
		    		log.error("Unable to perform assembly. Terminating early.");
	    		return -1;
	    	}	    	
	    	String instructions = getRealignmentScript(allEvidence);
	    	if (instructions != null && instructions.length() > 0) {
	    		log.error("Please realign intermediate fastq files. Suggested command-line for alignment is:\n" +
	    				"##################################\n"+
	    				instructions +
	    				"##################################");
		    	log.error("Please rerun after alignments have been performed.");
		    	if (SCRIPT != null) {
		    		FileWriter writer = null;
		    		try {
		    			writer = new FileWriter(SCRIPT);
		    			writer.write(instructions);
		    		} finally {
		    			if (writer != null) writer.close(); 
		    		}
		    		log.error("Realignment script has been written to ", SCRIPT);
		    	}
		    	return -1;
	    	}
	    	VariantCaller caller = null;
	    	try {
	    		caller = new VariantCaller(getContext(), OUTPUT, allEvidence);
	    		caller.callBreakends();
	    		caller.annotateBreakpoints();
	    	} finally {
	    		if (caller != null) caller.close();
	    	}
    	} catch (IOException e) {
    		log.error(e);
    		throw new RuntimeException(e);
    	}
		return 0;
    }
	public static void main(String[] argv) {
        System.exit(new Idsv().instanceMain(argv));
    }
}
