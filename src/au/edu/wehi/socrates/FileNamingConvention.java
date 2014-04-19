package au.edu.wehi.socrates;

import java.io.File;
import java.io.IOException;

/**
 * Encapsulates file naming conventions of intermediate and output files.
 * @author cameron.d
 *
 */
public class FileNamingConvention {
	private static final String COMMON_INITIAL_SUFFIX = ".socrates";
	private static final String FORMAT_SV_BAM = "%s" + COMMON_INITIAL_SUFFIX + ".sv.bam";
	private static final String FORMAT_ISV_BAM_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s.sv.bam";
	private static final String FORMAT_MATE_BAM = "%s" + COMMON_INITIAL_SUFFIX + ".svmate.bam";
	private static final String FORMAT_MATE_BAM_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s.svmate.bam";
	private static final String FORMAT_METRICS = "%s" + COMMON_INITIAL_SUFFIX + ".insertsize.txt";
	private static String GetStem(File input) throws IOException {
		String full = input.getAbsoluteFile().toString();
		if (full.contains(COMMON_INITIAL_SUFFIX)) {
			full = full.substring(0, full.indexOf(COMMON_INITIAL_SUFFIX));
		}
		return full;
	}
	public static File GetSVBam(File input) throws IOException {
		return new File(String.format(FORMAT_SV_BAM, GetStem(input)));
	}
	public static File GetSVBamForChr(File input, String chromosome) throws IOException {
		return new File(String.format(FORMAT_ISV_BAM_PER_CHR, GetStem(input), chromosome));
	}
	public static File GetMateBam(File input) throws IOException {
		return new File(String.format(FORMAT_MATE_BAM, GetStem(input)));
	}
	public static File GetMateBamForChr(File input, String chromosome) throws IOException {
		return new File(String.format(FORMAT_MATE_BAM_PER_CHR, GetStem(input), chromosome));
	}
	public static File GetMetrics(File input) throws IOException {
		return new File(String.format(FORMAT_METRICS, GetStem(input)));
	}
	//public static File getWorkingDirectory(File input) throws IOException {
	//	return new File(GetStem(input)).getParentFile();
	//}
	private FileNamingConvention() {}
}
