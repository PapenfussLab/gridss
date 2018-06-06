package gridss.cmdline;

import java.io.File;
import java.io.FileNotFoundException;
import java.nio.file.Path;

import org.broadinstitute.barclay.argparser.Argument;

import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.picard.TwoBitBufferedReferenceSequenceFile;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgram;

public abstract class ReferenceCommandLineProgram extends CommandLineProgram {
	private static final Log log = Log.getInstance(ReferenceCommandLineProgram.class);
	// --- intermediate file parameters ---
    @Argument(doc = "Directory to place intermediate results directories. Default location is the same directory"
    		+ " as the associated input or output file.", optional = true)
    public File WORKING_DIR = null;
    @Argument(doc = "Ignore reads marked as duplicates.", optional = true)
    public boolean IGNORE_DUPLICATES = true;
    private FileSystemContext fsc;
	private ReferenceLookup reference;
	public ReferenceLookup getReference() {
		IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
		if (reference == null) {
			ensureSequenceDictionary(REFERENCE_SEQUENCE, this);
			try {
				reference = new TwoBitBufferedReferenceSequenceFile(new IndexedFastaSequenceFile(REFERENCE_SEQUENCE));
			} catch (FileNotFoundException e) {
				String msg = String.format("Missing reference genome %s", REFERENCE_SEQUENCE);
				log.error(msg);
				throw new RuntimeException(msg);
			}
		}
		if (reference.getSequenceDictionary() == null) {
			Path referenceFile = REFERENCE_SEQUENCE.toPath();
			log.info("Attempting to create sequence dictionary for " + REFERENCE_SEQUENCE);
			Path dictPath = referenceFile.resolveSibling(referenceFile.getFileName().toString() + htsjdk.samtools.util.IOUtil.DICT_FILE_EXTENSION);
			picard.sam.CreateSequenceDictionary csd = new picard.sam.CreateSequenceDictionary();
			csd.instanceMain(new String[] {
				"O=" + dictPath.toFile(),
				"R=" + REFERENCE_SEQUENCE.getPath()
			});
		}
		return reference;
	}
	/**
	 * Ensures that a sequence dictionary exists for the given reference
	 * @param referenceFile reference genome fasta
	 * @param invoking program
	 */
	public static void ensureSequenceDictionary(File referenceFile, CommandLineProgram program) {
		try {
			ReferenceSequenceFile rsf = new FastaSequenceFile(referenceFile, false);
			Path path = referenceFile.toPath().toAbsolutePath();
			if (rsf.getSequenceDictionary() == null) {
				log.info("Attempting to create sequence dictionary for " + referenceFile);
				Path dictPath = path.resolveSibling(path.getFileName().toString() + htsjdk.samtools.util.IOUtil.DICT_FILE_EXTENSION);
				picard.sam.CreateSequenceDictionary csd = new picard.sam.CreateSequenceDictionary();
				if (program != null) {
					CommandLineProgramHelper.copyInputs(program, csd);
				}
				csd.instanceMain(new String[] {
					"OUTPUT=" + dictPath.toFile(),
					"R=" + referenceFile.getPath()
				});
			}
			rsf.close();
		} catch (Exception e) {
			log.error("Sequence dictionary creation failed. Please create using picard CreateSequenceDictionary.", e);
		}
	}
	public void setReference(ReferenceLookup ref) {
		this.reference = ref;
	}
	public FileSystemContext getFileSystemContext() {
		if (fsc == null) {
			fsc = new FileSystemContext(TMP_DIR.get(0), WORKING_DIR, MAX_RECORDS_IN_RAM);
		}
		return fsc;
	}
	public void setFileSystemContext(FileSystemContext fsc) {
		this.fsc = fsc;
	}
	@Override
	protected String[] customCommandLineValidation() {
		String[] val = referenceCustomCommandLineValidation();
		if (val != null) return val;
		return super.customCommandLineValidation();
	}
	public String[] referenceCustomCommandLineValidation() {
		if (referenceRequired()) {
			if (REFERENCE_SEQUENCE == null) {
	            return new String[]{"Must have a non-null REFERENCE_SEQUENCE"};
	        }
		}
		return null;
	}
	public boolean referenceRequired() { return true; }
	/**
	 * Copies the command line inputs to the given program
	 * @param cmd program to set command line for
	 */
	public void copyInputs(CommandLineProgram cmd) {
		CommandLineProgramHelper.copyInputs(this, cmd);
		
		if (cmd instanceof ReferenceCommandLineProgram) {
			ReferenceCommandLineProgram prog = (ReferenceCommandLineProgram) cmd;
			prog.WORKING_DIR = WORKING_DIR;
			prog.IGNORE_DUPLICATES = IGNORE_DUPLICATES;
			prog.fsc = fsc;
			prog.reference = reference;
		}
	}
}
