package gridss.cmdline;

import au.edu.wehi.idsv.FileSystemContext;
import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.picard.TwoBitBufferedReferenceSequenceFile;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.barclay.argparser.Argument;
import picard.cmdline.CommandLineProgram;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Path;
import java.util.List;

public abstract class ReferenceCommandLineProgram extends CommandLineProgram {
	private static final Log log = Log.getInstance(ReferenceCommandLineProgram.class);
	public static final List<String> BWA_COMMAND_LINE = ImmutableList.of(
			"bwa",
			"mem",
			"-K", "10000000",
			"-L", "0,0",
			"-t", "%3$d",
			"%2$s",
			"%1$s");
	public static final List<String> BOWTIE2_COMMAND_LINE = ImmutableList.of(
			"bowtie2",
			"--threads", "%3$d",
			"--local",
			"--mm",
			"--reorder",
			"-x", "%2$s",
			"-U", "%1$s");
	// --- intermediate file parameters ---
    @Argument(doc = "Directory to place intermediate results directories. Default location is the same directory"
    		+ " as the associated input or output file.", optional = true)
    public File WORKING_DIR = null;
    @Argument(doc = "Ignore reads marked as duplicates.", optional = true)
    public boolean IGNORE_DUPLICATES = true;
    private FileSystemContext fsc;
	private ReferenceLookup reference;
	public ReferenceLookup getReference() {
		if (reference == null) {
			IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
			ensureSequenceDictionary(REFERENCE_SEQUENCE);
			try {
				reference = new TwoBitBufferedReferenceSequenceFile(new IndexedFastaSequenceFile(REFERENCE_SEQUENCE));
			} catch (FileNotFoundException e) {
				String msg = String.format("Missing reference genome %s", REFERENCE_SEQUENCE);
				log.error(msg);
				throw new RuntimeException(msg);
			}
		}
		return reference;
	}
	/**
	 * Ensures that a sequence dictionary exists for the given reference
	 * @param referenceFile reference genome fasta
	 * @return true if a sequence dictionary was found, false if it had to be created
	 */
	public static boolean ensureSequenceDictionary(File referenceFile) {
		try (ReferenceSequenceFile rsf = new FastaSequenceFile(referenceFile, false)) {
			Path path = referenceFile.toPath().toAbsolutePath();
			Path dictPath = path.resolveSibling(path.getFileName().toString() + FileExtensions.DICT);
			if (!dictPath.toFile().exists()) {
				log.info("Attempting to create sequence dictionary for " + referenceFile);
				CommandLineProgramHelper cmd = new CommandLineProgramHelper(new picard.sam.CreateSequenceDictionary());
				cmd.addArg("OUTPUT", dictPath.toFile());
				cmd.addArg("R", referenceFile.getPath());
				cmd.run();
				return true;
			}
		} catch (Exception e) {
			log.error("Sequence dictionary creation failed. Please create using picard CreateSequenceDictionary.", e);
		}
		return false;
	}
	public void setReference(ReferenceLookup ref) {
		this.reference = ref;
	}
	public void setReference(File ref) {
		this.REFERENCE_SEQUENCE = ref;
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
	public void ensureDictionaryMatches(File input) throws IOException {
		SAMSequenceDictionary dictionary = getReference().getSequenceDictionary();
		ensureDictionaryMatches(input, dictionary, REFERENCE_SEQUENCE);
	}
	public static void ensureDictionaryMatches(File input, SAMSequenceDictionary referenceDictionary, File reference) throws IOException {
		final SamReaderFactory samFactory = SamReaderFactory.makeDefault().referenceSequence(reference);
		SamReader reader = null;
		reader = samFactory.open(input);
		try {
			final SAMSequenceDictionary samDictionary = reader.getFileHeader().getSequenceDictionary();
			if (samDictionary == null || samDictionary.isEmpty()) {
				String message = String.format("Missing @SQ sequencing dictionary header lines in %s. "
						+ " Are you sure this is a SAM/BAM/CRAM file? If so, make sure the @SQ headers are correct.", input);
				log.error(message);
				throw new RuntimeException(message);
			}
			SequenceUtil.assertSequenceDictionariesEqual(samDictionary, referenceDictionary, input, reference);
		} catch (htsjdk.samtools.util.SequenceUtil.SequenceListsDifferException e) {
			log.error("Reference genome used by ", input, " does not match reference genome ", reference, ". ",
					"The reference supplied must match the reference used for every input.");
			throw e;
		} finally {
			if (reader != null) reader.close();
		}
	}
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
