package au.edu.wehi.socrates;

import java.io.File;
import java.io.IOException;

import au.edu.wehi.socrates.FileNamingConvention.IntermediateBamType;

/**
 * Encapsulates file naming conventions of intermediate and output files.
 * @author cameron.d
 *
 */
public class FileNamingConvention {
	private static final String FORMAT_INTERMEDIATE_BAM = "%s_%s.bam";
	private static final String FORMAT_INTERMEDIATE_BAM_PER_CHR = "%s_%s_%s.bam";
	private static final String FORMAT_METRICS = "%s_insertsize.txt";
	public enum IntermediateBamType {
		/**
		 * Soft clipped reads
		 */
		SC,
		/**
		 * Open-ended anchor read pairs. Only one read in each pair is mapped.
		 */
		OEA,
		/**
		 * Discordant read pairs
		 */
		DP,
		/**
		 * Mate-coordinate sorted unmapped OEA reads, and DP read pairs.
		 */
		MATE
	}
	private static String GetStem(File input) throws IOException {
		String full = input.getAbsoluteFile().toString();
		return full;
	}
	public static File GetIntermediateBam(File input, IntermediateBamType type) throws IOException {
		return new File(String.format(FORMAT_INTERMEDIATE_BAM, GetStem(input), type.toString().toLowerCase()));
	}
	public static File GetIntermediateBamForChr(File input, IntermediateBamType type, String chromosome) throws IOException {
		return new File(String.format(FORMAT_INTERMEDIATE_BAM_PER_CHR, GetStem(input), type.toString().toLowerCase(), chromosome));
	}
	public static File GetMetrics(File input) throws IOException {
		return new File(String.format(FORMAT_METRICS, GetStem(input)));
	}
	//public static File getWorkingDirectory(File input) throws IOException {
	//	return new File(GetStem(input)).getParentFile();
	//}
	private FileNamingConvention() {}
}
