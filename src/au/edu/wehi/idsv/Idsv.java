package au.edu.wehi.idsv;

import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.EnumSet;
import java.util.List;

import picard.cmdline.Option;
import picard.cmdline.Usage;

import com.google.common.collect.Lists;

/**
 * Extracts structural variation evidence and assembles breakends
 * @author Daniel Cameron
 *
 */
public class Idsv extends CommandLineProgram {
	private Log log = Log.getInstance(Idsv.class);
	private static final String PROGRAM_VERSION = "0.1";
    // The following attributes define the command-line arguments
    @Usage
    public String USAGE = getStandardUsagePreamble() + "Calls structural variations from NGS sequencing data. " + PROGRAM_VERSION;
    
    @Option(doc = "Processing steps to execute",
            optional = true)
	public EnumSet<ProcessStep> STEPS = ProcessStep.ALL_STEPS;
    @Override
	protected int doWork() {
    	try {
	    	ensureDictionariesMatch();

	    	List<SAMEvidenceSource> samEvidence = Lists.newArrayList();
	    	for (File f : INPUT) {
	    		SAMEvidenceSource sref = new SAMEvidenceSource(getContext(), f, false);
	    		if (!sref.isComplete(ProcessStep.CALCULATE_METRICS)
	    			|| !sref.isComplete(ProcessStep.EXTRACT_SOFT_CLIPS)
	    			|| !sref.isComplete(ProcessStep.EXTRACT_READ_PAIRS)
	    			|| !sref.isComplete(ProcessStep.EXTRACT_READ_MATES)
	    			|| !sref.isComplete(ProcessStep.SORT_READ_MATES)) {
	    			sref.completeSteps(STEPS);
	    		}
	    		samEvidence.add(sref);
	    	}
	    	for (File f : INPUT_TUMOUR) {
	    		SAMEvidenceSource sref = new SAMEvidenceSource(getContext(), f, true);
	    		if (!sref.isComplete(ProcessStep.CALCULATE_METRICS)
	    			|| !sref.isComplete(ProcessStep.EXTRACT_SOFT_CLIPS)
	    			|| !sref.isComplete(ProcessStep.EXTRACT_READ_PAIRS)
	    			|| !sref.isComplete(ProcessStep.EXTRACT_READ_MATES)
	    			|| !sref.isComplete(ProcessStep.SORT_READ_MATES)) {
	    			sref.completeSteps(STEPS);
	    		}
	    		samEvidence.add(sref);
	    	}
	    	log.debug("Evidence extraction complete.");

	    	AssemblyEvidenceSource aes = new AssemblyEvidenceSource(getContext(), samEvidence, OUTPUT);
	    	aes.ensureAssembled();
	    	List<EvidenceSource> allEvidence = Lists.newArrayList();
	    	allEvidence.add(aes);
	    	allEvidence.addAll(samEvidence);
	    	
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
