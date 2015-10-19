package au.edu.wehi.validation;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

/**
 * Validations deletion calls with long read support
 * @author Daniel Cameron
 *
 */
@CommandLineProgramProperties(
		usage = "Validations deletions 50bp or more from long read alignments.",
		usageShort = "Validations deletions 50bp or more from long read alignments.")
public class ValidateDeletions extends CommandLineProgram {
    @Option(doc="Deletions in BEDPE format", shortName=StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File BEDPE = null;
    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;
    @Option(doc="Long read BAM file", shortName="LR")
    public List<File> LONG_READS;
    //@Option(shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="Reference used for alignment")
    //public File REFERENCE;
    @Option(doc="Minimum left of a long read soft clip before it is considered supporting", shortName="MSC", optional=true)
    public int MIN_SOFT_CLIP_LENGTH = 30;
    @Option(doc="Error margin allowed between the soft clip position and the call position", shortName="SCM", optional=true)
    public int SOFT_CLIP_MARGIN = 30;
    @Option(doc="Number of bases to expand the called deletion when counting spanning deletion events", shortName="SWS", optional=true)
    public int SPANNING_WINDOW_SIZE = 50;
    @Option(doc="Minimum size of a deletion before it is considered to contribute to a spanning event."
    		+ " This is used to remove PacBio short indel from the deletion signal.", shortName="MDS", optional=true)
    public int MIN_DELETION_LENGTH = 4;
    private static final String EOF_INDICATOR = "EOF";
    protected int doWork() {
        try {
        	//IOUtil.assertFileIsReadable(REFERENCE);
        	
        	List<ArrayBlockingQueue<String>> queues = new ArrayList<ArrayBlockingQueue<String>>();
        	for (int i = 0; i < LONG_READS.size() + 1; i++) {
        		queues.add(new ArrayBlockingQueue<String>(1024));
        	}
        	new Thread(() -> fill(BEDPE, queues.get(0))).start();
        	for (int i = 0; i < LONG_READS.size(); i++) {
        		LongReadSupportFinder finder = new LongReadSupportFinder(LONG_READS.get(i));
        		ArrayBlockingQueue<String> input = queues.get(i);
        		ArrayBlockingQueue<String> output = queues.get(i + 1);
        		new Thread(() -> annotate(finder, input, output)).start();
        	}
        	ArrayBlockingQueue<String> output = queues.get(queues.size() - 1);
        	BufferedWriter writer = new BufferedWriter(new FileWriter(OUTPUT));
        	while (true) {
        		String line = output.take();
        		if (line == EOF_INDICATOR) break;
        		writer.write(line);
        		writer.write('\n');
        	}
        	writer.close();
        } catch (Exception e) {
			e.printStackTrace();
			return 1;
		}
        return 0;
    }
    public static void fill(File bedpe, ArrayBlockingQueue<String> output) {
    	try {
	    	for (String line : Files.readAllLines(bedpe.toPath())) {
	    		output.put(line);
	    	}
	    	output.put(EOF_INDICATOR);
    	} catch (Exception e) {
			e.printStackTrace();
		}
    }
    public void annotate(LongReadSupportFinder finder, ArrayBlockingQueue<String> in, ArrayBlockingQueue<String> out) {
    	try {
    		while (true) {
	    		String line = in.take();
	    		if (line == EOF_INDICATOR) break;
	    		StringBuilder sb = new StringBuilder(line);
	    		if (line.startsWith("#")) {
	    			sb.append(BedpeDeletion.FS);
	    			sb.append("SpanningDeletionSize");
	    			sb.append(BedpeDeletion.FS);
	    			sb.append("StartClipLocation");
	    			sb.append(BedpeDeletion.FS);
	    			sb.append("EndClipLocation");
	    		} else {
		    		BedpeDeletion del = new BedpeDeletion(line);
		    		LongReadSupportLevel support = finder.evaluateDeletion(
		    				del.chrom1, del.start1, del.end1, del.start2, del.end2,
		    				SOFT_CLIP_MARGIN, MIN_SOFT_CLIP_LENGTH, SPANNING_WINDOW_SIZE, MIN_DELETION_LENGTH);
		    		if (support == null) {
		    			sb.append(BedpeDeletion.FS);
		    			sb.append("ERROR");
		        		sb.append(BedpeDeletion.FS);
		        		sb.append("ERROR");
		        		sb.append(BedpeDeletion.FS);
		        		sb.append("ERROR");
		    		} else {
		    			sb.append(BedpeDeletion.FS);
		    			sb.append(support.spanningAlignments.toString());
		    			sb.append(BedpeDeletion.FS);
		        		sb.append(support.startClipLocations.toString());
		        		sb.append(BedpeDeletion.FS);
		        		sb.append(support.endClipLocations.toString());
		    		}
	    		}
	    		out.put(sb.toString());
	    	}
	    	out.put(EOF_INDICATOR);
    	} catch (Exception e) {
			e.printStackTrace();
		}
    }
	public static void main(String[] argv) {
	    System.exit(new ValidateDeletions().instanceMain(argv));
	}
}
