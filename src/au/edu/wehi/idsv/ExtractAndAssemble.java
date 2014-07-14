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
	    	List<SAMEvidenceSource> evidence = Lists.newArrayList();
	    	for (File f : INPUT) {
	    		SAMEvidenceSource sref = new SAMEvidenceSource(getContext(), f, false);
	    		sref.ensureEvidenceExtracted();
	    		evidence.add(sref);
	    	}
	    	for (File f : INPUT_TUMOUR) {
	    		SAMEvidenceSource sref = new SAMEvidenceSource(getContext(), f, true);
	    		sref.ensureEvidenceExtracted();
	    		evidence.add(sref);
	    	}
	    	log.info("Evidence extraction complete.");
	    	log.info("Starting assembly");
	    	AssemblyEvidenceSource aes = new AssemblyEvidenceSource(getContext(), evidence, OUTPUT);
	    	aes.ensureAssembled();
	    	log.info("Assembly complete");
	    	printAssemblyInstructions(Iterables.concat(evidence, ImmutableList.of(aes)));
    	} catch (IOException e) {
    		log.error(e);
    		throw new RuntimeException(e);
    	}
		return 0;
    }
    private void printAssemblyInstructions(Iterable<EvidenceSource> it) {
    	int realignmentCount = 0;
    	StringBuilder sb = new StringBuilder("Please realign intermediate fastq files. Suggested command-line for alignment is:\n");
    	sb.append("########## Start Recommended Aligner Commands ##########\n");
    	for (EvidenceSource source : it) {
    		if (!source.isRealignmentComplete()) {
    			realignmentCount++;
    			sb.append(source.getRealignmentScript());
    		}
    	}
    	sb.append("########## End Recommended Aligner Commands ##########\n");
    	if (realignmentCount > 0) {
    		log.info(sb);
    	}
    }
	public static void main(String[] argv) {
        System.exit(new ExtractAndAssemble().instanceMain(argv));
    }
}
