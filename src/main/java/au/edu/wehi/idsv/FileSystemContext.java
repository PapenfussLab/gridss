package au.edu.wehi.idsv;

import java.io.File;
import java.io.IOException;

public class FileSystemContext {
	private final File tempDir;
	private final File workingDir;
	private final int maxRecordsInRam;
	public FileSystemContext(
			File tempDir,
			File workingDir,
			int maxRecordsInRam) {
		if (tempDir.exists() && !tempDir.isDirectory()) throw new IllegalArgumentException(String.format("temp directory %s is not a directory", tempDir));
		if (workingDir != null && workingDir.exists() && !workingDir.isDirectory()) throw new IllegalArgumentException(String.format("working directory %s is not a directory", workingDir));
		this.tempDir = tempDir;
		this.workingDir = workingDir;
		this.maxRecordsInRam = maxRecordsInRam;
	}
	public FileSystemContext(
			File tempDir,
			int maxRecordsInRam) {
		this(tempDir, null, maxRecordsInRam);
	}
	public static File getWorkingFileFor(File file) {
		return getWorkingFileFor(file, "tmp.");
	}
	public static File getWorkingFileFor(File file, String workingPrefix) {
		return new File(file.getParent(), workingPrefix + file.getName());
	}
	public File getTemporaryDirectory() {
		return tempDir;
	}
	/**
	 * Maximum number of buffered records per file input or output stream.
	 * @see Corresponds to {@link picard.cmdline.CommandLineProgram.MAX_RECORDS_IN_RAM}
	 * @return maximum records in memory
	 */
	public int getMaxBufferedRecordsPerFile() {
		return maxRecordsInRam;
	}
	private static final String SAM_SUFFIX = ".bam";
	private static final String VCF_SUFFIX = ".vcf";
	private static final String COMMON_INITIAL_SUFFIX = ".idsv";
	private static final String INTERMEDIATE_DIR_SUFFIX = COMMON_INITIAL_SUFFIX + ".working";
	private static final String FORMAT_SC_SAM = "%s" + COMMON_INITIAL_SUFFIX + ".sc" + SAM_SUFFIX;
	private static final String FORMAT_SC_SAM_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s.sc" + SAM_SUFFIX;
	private static final String FORMAT_RP_SAM = "%s" + COMMON_INITIAL_SUFFIX + ".rp" + SAM_SUFFIX;
	private static final String FORMAT_RP_SAM_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s.rp" + SAM_SUFFIX;
	private static final String FORMAT_MATE_SAM = "%s" + COMMON_INITIAL_SUFFIX + ".rpmate" + SAM_SUFFIX;
	private static final String FORMAT_MATE_SAM_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s.rpmate" + SAM_SUFFIX;
	private static final String FORMAT_MATE_SAM_UNSORTED = "%s" + COMMON_INITIAL_SUFFIX + ".rpmate.unsorted" + SAM_SUFFIX;
	private static final String FORMAT_MATE_SAM_UNSORTED_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s.rpmate.unsorted" + SAM_SUFFIX;
	private static final String FORMAT_INSERT_SIZE_METRICS = "%s" + COMMON_INITIAL_SUFFIX + ".metrics.insertsize.txt";
	private static final String FORMAT_IDSV_METRICS = "%s" + COMMON_INITIAL_SUFFIX + ".metrics.idsv.txt";
	private static final String FORMAT_REALIGN_FASTQ = "%s" + COMMON_INITIAL_SUFFIX + ".realign.fq";
	private static final String FORMAT_REALIGN_FASTQ_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s.realign.fq";
	private static final String FORMAT_REALIGN_SAM = "%s" + COMMON_INITIAL_SUFFIX + ".realign" + SAM_SUFFIX;
	private static final String FORMAT_REALIGN_SAM_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s.realign" + SAM_SUFFIX;
	private static final String FORMAT_ASSEMBLY_RAW = "%s" + COMMON_INITIAL_SUFFIX + ".breakend" + SAM_SUFFIX;
	private static final String FORMAT_ASSEMBLY_RAW_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s.breakend" + SAM_SUFFIX;
	private static final String FORMAT_BREAKPOINT_VCF_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s-%s.breakpoint" + VCF_SUFFIX;
	private static final String FORMAT_BREAKPOINT_VCF = "%s" + COMMON_INITIAL_SUFFIX + ".breakpoint" + VCF_SUFFIX;
	private static final String FORMAT_FLAG_FILE = "%s" + COMMON_INITIAL_SUFFIX + ".%s";
	private static final String FORMAT_SC_REMOTE_SAM = "%s" + COMMON_INITIAL_SUFFIX + ".scremote" + SAM_SUFFIX;
	private static final String FORMAT_SC_REMOTE_SAM_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s.scremote" + SAM_SUFFIX;
	private static final String FORMAT_SC_REMOTE_SAM_UNSORTED = "%s" + COMMON_INITIAL_SUFFIX + ".scremote.unsorted" + SAM_SUFFIX;
	private static final String FORMAT_SC_REMOTE_SAM_UNSORTED_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s.scremote.unsorted" + SAM_SUFFIX;;
	private static final String FORMAT_REALIGN_REMOTE_SAM_UNSORTED = "%s" + COMMON_INITIAL_SUFFIX + ".realignremote.unsorted" + SAM_SUFFIX;
	private static final String FORMAT_REALIGN_REMOTE_SAM_PER_CHR_UNSORTED = "%s" + COMMON_INITIAL_SUFFIX + ".%s.realignremote.unsorted" + SAM_SUFFIX;
	private static final String FORMAT_REALIGN_REMOTE_SAM = "%s" + COMMON_INITIAL_SUFFIX + ".realignremote" + SAM_SUFFIX;
	private static final String FORMAT_REALIGN_REMOTE_SAM_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s.realignremote" + SAM_SUFFIX;
	private static final String FORMAT_ASSEMBLY = "%s" + COMMON_INITIAL_SUFFIX + ".assembly" + SAM_SUFFIX;
	private static final String FORMAT_ASSEMBLY_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s.assembly" + SAM_SUFFIX;
	private static final String FORMAT_ASSEMBLY_MATE = "%s" + COMMON_INITIAL_SUFFIX + ".assemblymate" + SAM_SUFFIX;
	private static final String FORMAT_ASSEMBLY_MATE_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s.assemblymate" + SAM_SUFFIX;
	
	/**
	 * Gets the idsv intermediate working directory for the given input
	 * @param input
	 * @return
	 */
	public File getWorkingDirectory(File file) {
		if (workingDir != null) {
			return workingDir;
		}
		return getSource(file).getParentFile();
	}
	/**
	 * Ensures that the intermediate directory for the given file exists
	 * @param file file
	 * @return intermediate directory
	 */
	public File getIntermediateDirectory(File file) {
		File input = getSource(file);
		File workingDir = new File(getWorkingDirectory(file), input.getName() + INTERMEDIATE_DIR_SUFFIX);
		if (!workingDir.exists()) {
			workingDir.mkdir();
		}
		return workingDir;
	}
	private static File getSource(File file) {
		String source = file.getAbsolutePath();
		if (source.contains(INTERMEDIATE_DIR_SUFFIX)) {
			source = source.substring(0, source.indexOf(INTERMEDIATE_DIR_SUFFIX));
		}
		File result = new File(source); 
		return result;
	}
	/**
	 * Gets the working directory stem for the given input
	 * @param file input bam or idsv intermediate file
	 * @return working location filename prefix for the input file 
	 * @throws IOException
	 */
	private String getStem(File file) {
		File stem = new File(getIntermediateDirectory(file), getSource(file).getName());
		return stem.getAbsolutePath();

	}
	public File getReadPairBam(File input) {
		return new File(String.format(FORMAT_RP_SAM, getStem(input)));
	}
	public File getReadPairBamForChr(File input, String chromosome) {
		return new File(String.format(FORMAT_RP_SAM_PER_CHR, getStem(input), chromosome));
	}
	public File getSoftClipBam(File input) {
		return new File(String.format(FORMAT_SC_SAM, getStem(input)));
	}
	public File getSoftClipBamForChr(File input, String chromosome) {
		return new File(String.format(FORMAT_SC_SAM_PER_CHR, getStem(input), chromosome));
	}
	public File getMateBam(File input) {
		return new File(String.format(FORMAT_MATE_SAM, getStem(input)));
	}
	public File getMateBamForChr(File input, String chromosome) {
		return new File(String.format(FORMAT_MATE_SAM_PER_CHR, getStem(input), chromosome));
	}
	public File getMateBamUnsorted(File input) {
		return new File(String.format(FORMAT_MATE_SAM_UNSORTED, getStem(input)));
	}
	public File getMateBamUnsortedForChr(File input, String chromosome) {
		return new File(String.format(FORMAT_MATE_SAM_UNSORTED_PER_CHR, getStem(input), chromosome));
	}
	public File getAssemblyRawBam(File input) {
		return new File(String.format(FORMAT_ASSEMBLY_RAW, getStem(input)));
	}
	public File getAssemblyRawBamForChr(File input, String chromosome) {
		return new File(String.format(FORMAT_ASSEMBLY_RAW_PER_CHR, getStem(input), chromosome));
	}
	public File getAssembly(File input) {
		return new File(String.format(FORMAT_ASSEMBLY, getStem(input)));
	}
	public File getAssemblyForChr(File input, String chromosome) {
		return new File(String.format(FORMAT_ASSEMBLY_PER_CHR, getStem(input), chromosome));
	}
	public File getAssemblyMate(File input) {
		return new File(String.format(FORMAT_ASSEMBLY_MATE, getStem(input)));
	}
	public File getAssemblyMateForChr(File input, String chromosome) {
		return new File(String.format(FORMAT_ASSEMBLY_MATE_PER_CHR, getStem(input), chromosome));
	}
	public File getBreakpointVcf(File input) {
		return new File(String.format(FORMAT_BREAKPOINT_VCF, getStem(input)));
	}
	public File getBreakpointVcf(File input, String chromosome1, String chromosome2) {
		return new File(String.format(FORMAT_BREAKPOINT_VCF_PER_CHR, getStem(input), chromosome1, chromosome2));
	}
	public File getInsertSizeMetrics(File input) {
		return new File(String.format(FORMAT_INSERT_SIZE_METRICS, getStem(input)));
	}
	public File getIdsvMetrics(File input) {
		return new File(String.format(FORMAT_IDSV_METRICS, getStem(input)));
	}
	public File getRealignmentBam(File input) {
		return new File(String.format(FORMAT_REALIGN_SAM, getStem(input)));
	}
	public File getRealignmentBamForChr(File input, String chromosome) {
		return new File(String.format(FORMAT_REALIGN_SAM_PER_CHR, getStem(input), chromosome));
	}
	public File getRealignmentFastq(File input) {
		return new File(String.format(FORMAT_REALIGN_FASTQ, getStem(input)));
	}
	public File getRealignmentFastqForChr(File input, String chromosome) {
		return new File(String.format(FORMAT_REALIGN_FASTQ_CHR, getStem(input), chromosome));
	}
	public File getFlagFile(File input, String flagname) {
		return new File(String.format(FORMAT_FLAG_FILE, getStem(input), flagname));
	}
	// Remote sorted
	public File getSoftClipRemoteBam(File input) {
		return new File(String.format(FORMAT_SC_REMOTE_SAM, getStem(input)));
	}
	public File getSoftClipRemoteBamForChr(File input, String chromosome) {
		return new File(String.format(FORMAT_SC_REMOTE_SAM_PER_CHR, getStem(input), chromosome));
	}
	public File getSoftClipRemoteUnsortedBam(File input) {
		return new File(String.format(FORMAT_SC_REMOTE_SAM_UNSORTED, getStem(input)));
	}
	public File getSoftClipRemoteUnsortedBamForChr(File input, String chromosome) {
		return new File(String.format(FORMAT_SC_REMOTE_SAM_UNSORTED_PER_CHR, getStem(input), chromosome));
	}
	public File getRealignmentRemoteUnsortedBam(File input) {
		return new File(String.format(FORMAT_REALIGN_REMOTE_SAM_UNSORTED, getStem(input)));
	}
	public File getRealignmentRemoteUnsortedBamForChr(File input, String chromosome) {
		return new File(String.format(FORMAT_REALIGN_REMOTE_SAM_PER_CHR_UNSORTED, getStem(input), chromosome));
	}
	public File getRealignmentRemoteBam(File input) {
		return new File(String.format(FORMAT_REALIGN_REMOTE_SAM, getStem(input)));
	}
	public File getRealignmentRemoteBamForChr(File input, String chromosome) {
		return new File(String.format(FORMAT_REALIGN_REMOTE_SAM_PER_CHR, getStem(input), chromosome));
	}
}
