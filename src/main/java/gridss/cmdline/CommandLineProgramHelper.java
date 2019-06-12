package gridss.cmdline;

import picard.cmdline.CommandLineProgram;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

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
		to.referenceSequence = from.referenceSequence;
		to.TMP_DIR = from.TMP_DIR;
		to.VALIDATION_STRINGENCY = from.VALIDATION_STRINGENCY;
		to.VERBOSITY = from.VERBOSITY;
	}
	public static List<String> getCommonArgs(CommandLineProgram program) {
		List<String> args = new ArrayList<>();
		args.add("COMPRESSION_LEVEL=" + Integer.toString(program.COMPRESSION_LEVEL));
		args.add("CREATE_INDEX=" + (program.CREATE_INDEX == null ? "null" : program.CREATE_INDEX.toString()));
		args.add("CREATE_MD5_FILE=" + Boolean.toString(program.CREATE_MD5_FILE));
		args.add("GA4GH_CLIENT_SECRETS=" + (program.GA4GH_CLIENT_SECRETS == null ? "null" : program.GA4GH_CLIENT_SECRETS));
		args.add("MAX_RECORDS_IN_RAM=" + (program.MAX_RECORDS_IN_RAM == null ? "null" : program.MAX_RECORDS_IN_RAM.toString()));
		args.add("QUIET=" + (program.QUIET == null ? "null" : Boolean.toString(program.QUIET)));
		if (program.referenceSequence.getReferenceFile() != null) {
			args.add("REFERENCE_SEQUENCE=" + program.referenceSequence.getReferenceFile().getPath());	
		}
		args.add("VALIDATION_STRINGENCY=" + (program.VALIDATION_STRINGENCY == null ? "null" : program.VALIDATION_STRINGENCY.toString()));
		args.add("VERBOSITY=" + (program.VERBOSITY == null ? "null" : program.VERBOSITY.toString()));
		if (program.TMP_DIR == null || program.TMP_DIR.size() == 0) {
			args.add("TMP_DIR=null");
		} else {
			for (File tmp : program.TMP_DIR) {
				args.add("TMP_DIR=" + tmp.getPath());
			}
		}
		return args;
	}
}
