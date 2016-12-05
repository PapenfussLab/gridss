package gridss.cmdline;

import java.io.FileNotFoundException;

import au.edu.wehi.idsv.picard.ReferenceLookup;
import au.edu.wehi.idsv.picard.TwoBitBufferedReferenceSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgram;

public abstract class ReferenceCommandLineProgram extends CommandLineProgram {
	private static final Log log = Log.getInstance(ReferenceCommandLineProgram.class);
	private ReferenceLookup reference;
	public ReferenceLookup getReference() {
		IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
		if (reference == null) {
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
	public void setReference(ReferenceLookup ref) {
		this.reference =  ref;
	}
	@Override
	protected String[] customCommandLineValidation() {
		String[] val = referenceCustomCommandLineValidation();
		if (val != null) return val;
		return super.customCommandLineValidation();
	}
	public String[] referenceCustomCommandLineValidation() {
		if (REFERENCE_SEQUENCE == null) {
            return new String[]{"Must have a non-null REFERENCE_SEQUENCE"};
        }
		return null;
	}
	/**
	 * Copies the command line inputs to the given program
	 * @param cmd program to set command line for
	 */
	public void copyInputs(CommandLineProgram cmd) {
		cmd.COMPRESSION_LEVEL = COMPRESSION_LEVEL;
		cmd.CREATE_INDEX = CREATE_INDEX;
		cmd.CREATE_MD5_FILE = CREATE_MD5_FILE;
		cmd.GA4GH_CLIENT_SECRETS = GA4GH_CLIENT_SECRETS;
		cmd.MAX_RECORDS_IN_RAM = MAX_RECORDS_IN_RAM;
		cmd.QUIET = QUIET;
		cmd.REFERENCE_SEQUENCE = REFERENCE_SEQUENCE;
		cmd.TMP_DIR = TMP_DIR;
		cmd.VALIDATION_STRINGENCY = VALIDATION_STRINGENCY;
		cmd.VERBOSITY = VERBOSITY;
	}
}
