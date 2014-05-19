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
	private static final String FORMAT_SV_BAM_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s.sv.bam";
	private static final String FORMAT_MATE_BAM = "%s" + COMMON_INITIAL_SUFFIX + ".svmate.bam";
	private static final String FORMAT_MATE_BAM_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s.svmate.bam";
	private static final String FORMAT_METRICS = "%s" + COMMON_INITIAL_SUFFIX + ".metrics.txt";
	private static final String FORMAT_REALIGN_BAM = "%s" + COMMON_INITIAL_SUFFIX + ".realign.bam";
	private static final String FORMAT_REALIGN_BAM_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s.realign.bam";
	private static final String FORMAT_BREAKEND_VCF_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s.breakend.vcf";
	private static final String FORMAT_BREAKEND_VCF = "%s" + COMMON_INITIAL_SUFFIX + ".breakend.vcf";
	private static final String FORMAT_BREAKPOINT_VCF_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s-%s.breakpoint.vcf";
	private static final String FORMAT_BREAKPOINT_VCF = "%s" + COMMON_INITIAL_SUFFIX + ".breakpoint.vcf";
	private static final String FORMAT_OUTPUT_VCF = "%s" + COMMON_INITIAL_SUFFIX + ".vcf";
	private static String GetStem(File input) throws IOException {
		String full = input.getAbsoluteFile().toString();
		if (full.contains(COMMON_INITIAL_SUFFIX)) {
			full = full.substring(0, full.indexOf(COMMON_INITIAL_SUFFIX));
		}
		return full;
	}
	public static File getSVBam(File input) throws IOException {
		return new File(String.format(FORMAT_SV_BAM, GetStem(input)));
	}
	public static File getSVBamForChr(File input, String chromosome) throws IOException {
		return new File(String.format(FORMAT_SV_BAM_PER_CHR, GetStem(input), chromosome));
	}
	public static File getMateBam(File input) throws IOException {
		return new File(String.format(FORMAT_MATE_BAM, GetStem(input)));
	}
	public static File getMateBamForChr(File input, String chromosome) throws IOException {
		return new File(String.format(FORMAT_MATE_BAM_PER_CHR, GetStem(input), chromosome));
	}
	public static File getBreakendVcf(File input) throws IOException {
		return new File(String.format(FORMAT_BREAKEND_VCF, GetStem(input)));
	}
	public static File getBreakendVcfForChr(File input, String chromosome) throws IOException {
		return new File(String.format(FORMAT_BREAKEND_VCF_PER_CHR, GetStem(input), chromosome));
	}
	public static File getRawCallVcf(File input) throws IOException {
		return new File(String.format(FORMAT_BREAKPOINT_VCF, GetStem(input)));
	}
	public static File getRawCallVcf(File input, String chromosome1, String chromosome2) throws IOException {
		return new File(String.format(FORMAT_BREAKPOINT_VCF_PER_CHR, GetStem(input), chromosome1, chromosome2));
	}
	public static File getOutputVcf(File input) throws IOException {
		return new File(String.format(FORMAT_OUTPUT_VCF, GetStem(input)));
	}
	public static File getMetrics(File input) throws IOException {
		return new File(String.format(FORMAT_METRICS, GetStem(input)));
	}
	public static File getRealignmentBam(File input) throws IOException {
		return new File(String.format(FORMAT_REALIGN_BAM, GetStem(input)));
	}
	public static File getRealignmentBamForChr(File input, String chromosome) throws IOException {
		return new File(String.format(FORMAT_REALIGN_BAM_PER_CHR, GetStem(input), chromosome));
	}
	//public static File getWorkingDirectory(File input) throws IOException {
	//	return new File(GetStem(input)).getParentFile();
	//}
	private FileNamingConvention() {}
}
