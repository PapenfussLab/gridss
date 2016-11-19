package gridss;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Locale;

import org.apache.commons.lang.NotImplementedException;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMEvidenceSource;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(
        usage = "Extracts reads and read pairs supporting putative structural variations.",
        usageShort = "Extracts reads and read pairs supporting putative structural variations."
)
public class AssembleBreakends extends CommandLineProgram {
	private static final Log log = Log.getInstance(AssembleBreakends.class);
    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input file", optional=false)
    public List<File> INPUT;
    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output file containing subset of input", optional=false)
    public File OUTPUT;
    @Override
	protected int doWork() {
		log.debug("Setting language-neutral locale");
    	java.util.Locale.setDefault(Locale.ROOT);
    	validateParameters();
    	try {
	    	ProcessingContext pc = getContext();
	    	List<SAMEvidenceSource> sources = asEvidenceSources();
	    	AssemblyEvidenceSource assembler = new AssemblyEvidenceSource(pc, sources, OUTPUT);
	    	assembler.assembleBreakends();
	    } catch (IOException e) {
			log.error(e);
			return -1;
		}
    	return 0;
    }
    private List<SAMEvidenceSource> asEvidenceSources() {
    	throw new NotImplementedException();
    }
    private ProcessingContext getContext() {
    	throw new NotImplementedException();
    }
	private void validateParameters() {
		for (File f : INPUT) {
			IOUtil.assertFileIsReadable(f);
		}
    	IOUtil.assertFileIsWritable(OUTPUT);
	}
	public static void main(String[] argv) {
        System.exit(new AssembleBreakends().instanceMain(argv));
    }
}
