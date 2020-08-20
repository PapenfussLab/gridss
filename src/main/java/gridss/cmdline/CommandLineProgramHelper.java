package gridss.cmdline;

import com.google.common.collect.Lists;
import picard.cmdline.CommandLineProgram;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class CommandLineProgramHelper {
	private final List<String> argName = new ArrayList<>();
	private final List<Object> argValue = new ArrayList<>();
	private final List<String> positionalArgs = new ArrayList<>();
	private final CommandLineProgram program;

	public CommandLineProgramHelper(CommandLineProgram program) {
		this.program = program;
	}
	public CommandLineProgram getProgram() {
		return program;
	}
	public void addArg(String name, Object value) {
		argName.add(name);
		argValue.add(value);
	}
	public int run() {
		return run(program, argName, argValue, positionalArgs);
	}
	public static int run(CommandLineProgram program, String[] args, String... positionalArgs) {
		List<String> argName = new ArrayList<>();
		List<Object> argValue = new ArrayList<>();
		for (int i = 0; i < args.length; i += 2) {
			argName.add(args[i]);
			argValue.add(i + 1 <= args.length ? null : args[i + 1]);
		}
		return run(program, argName, argValue, Lists.newArrayList(positionalArgs));
	}
	public static int run(CommandLineProgram program, List<String> argName, List<Object> argValue, List<String> positionalArgs) {
		if (argName.size() != argValue.size()) {
			throw new IllegalArgumentException("argName and argValue are different lengths");
		}
		boolean useLegacySyntax = Boolean.parseBoolean(System.getProperty("picard.useLegacyParser", "true"));
		List<String> args = new ArrayList<>();
		for (int i = 0 ; i < argName.size(); i++) {
			String name = argName.get(i);
			Object objValue = argValue.get(i);
			String strValue = objValue == null ? "null" : objValue.toString();
			if (useLegacySyntax) {
				args.add(name + "=" + strValue);
			} else {
				args.add("-" + name);
				args.add(strValue);
			}
		}
		args.addAll(positionalArgs);
		return program.instanceMain(args.toArray(new String[0]));
	}
	public void copyInputs(CommandLineProgram from) {
		copyInputs(from, program);
	}
	public void setCommonArgs(CommandLineProgram from) {
		addArg("COMPRESSION_LEVEL", program.COMPRESSION_LEVEL);
		addArg("CREATE_INDEX", (program.CREATE_INDEX == null ? "null" : program.CREATE_INDEX.toString()));
		addArg("CREATE_MD5_FILE", program.CREATE_MD5_FILE);
		addArg("GA4GH_CLIENT_SECRETS", (program.GA4GH_CLIENT_SECRETS == null ? "null" : program.GA4GH_CLIENT_SECRETS));
		addArg("MAX_RECORDS_IN_RAM", (program.MAX_RECORDS_IN_RAM == null ? "null" : program.MAX_RECORDS_IN_RAM.toString()));
		addArg("QUIET", (program.QUIET == null ? "null" : Boolean.toString(program.QUIET)));
		if (program.referenceSequence.getReferenceFile() != null) {
			addArg("REFERENCE_SEQUENCE", program.referenceSequence.getReferenceFile().getPath());
		}
		addArg("VALIDATION_STRINGENCY", (program.VALIDATION_STRINGENCY == null ? "null" : program.VALIDATION_STRINGENCY.toString()));
		addArg("VERBOSITY", (program.VERBOSITY == null ? "null" : program.VERBOSITY.toString()));
		if (program.TMP_DIR == null || program.TMP_DIR.size() == 0) {
			addArg("TMP_DIR", "null");
		} else {
			for (File tmp : program.TMP_DIR) {
				addArg("TMP_DIR", tmp.getPath());
			}
		}
	}
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
}
