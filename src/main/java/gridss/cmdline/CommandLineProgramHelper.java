package gridss.cmdline;

import picard.cmdline.CommandLineProgram;

public class CommandLineProgramHelper {
	/**
	 * Copies the command line inputs to the given program
	 * @param to program to set command line for
	 */
	public static void copyInputs(CommandLineProgram from, CommandLineProgram to) {
		to.COMPRESSION_LEVEL = from.COMPRESSION_LEVEL;
		to.CREATE_INDEX = from.CREATE_INDEX;
		to.CREATE_MD5_FILE = from.CREATE_MD5_FILE;
		to.GA4GH_CLIENT_SECRETS = from.GA4GH_CLIENT_SECRETS;
		to.MAX_RECORDS_IN_RAM = from.MAX_RECORDS_IN_RAM;
		to.QUIET = from.QUIET;
		to.REFERENCE_SEQUENCE = from.REFERENCE_SEQUENCE;
		to.TMP_DIR = from.TMP_DIR;
		to.VALIDATION_STRINGENCY = from.VALIDATION_STRINGENCY;
		to.VERBOSITY = from.VERBOSITY;
	}
}
