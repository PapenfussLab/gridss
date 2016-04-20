package au.edu.wehi.idsv;

import java.io.File;
import java.util.HashSet;
import java.util.Set;

public class FileSystemContext {
	private final File tempDir;
	private final File workingDir;
	private final int maxRecordsInRam;
	/**
	 * List of directories already created/known to exist.
	 * Used as a performance optimisation to reduce expensive samba I/O operations
	 */
	private Set<File> createdDirectories = new HashSet<File>();
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
	private static final String FORMAT_SC_SAM = "%1$s/%2$s" + COMMON_INITIAL_SUFFIX + ".sc" + SAM_SUFFIX;
	private static final String FORMAT_SC_SAM_PER_CHR = "%1$s/%3$s/%2$s" + COMMON_INITIAL_SUFFIX + ".%3$s.sc" + SAM_SUFFIX;
	private static final String FORMAT_RP_SAM = "%1$s/%2$s" + COMMON_INITIAL_SUFFIX + ".rp" + SAM_SUFFIX;
	private static final String FORMAT_RP_SAM_PER_CHR = "%1$s/%3$s/%2$s" + COMMON_INITIAL_SUFFIX + ".%3$s.rp" + SAM_SUFFIX;
	private static final String FORMAT_MATE_SAM = "%1$s/%2$s" + COMMON_INITIAL_SUFFIX + ".rpmate" + SAM_SUFFIX;
	private static final String FORMAT_MATE_SAM_PER_CHR = "%1$s/%3$s/%2$s" + COMMON_INITIAL_SUFFIX + ".%3$s.rpmate" + SAM_SUFFIX;
	private static final String FORMAT_MATE_SAM_UNSORTED = "%1$s/%2$s" + COMMON_INITIAL_SUFFIX + ".rpmate.unsorted" + SAM_SUFFIX;
	private static final String FORMAT_MATE_SAM_UNSORTED_PER_CHR = "%1$s/%3$s/%2$s" + COMMON_INITIAL_SUFFIX + ".%3$s.rpmate.unsorted" + SAM_SUFFIX;
	private static final String FORMAT_INSERT_SIZE_METRICS = "%1$s/%2$s" + COMMON_INITIAL_SUFFIX + ".metrics.insertsize.txt";
	private static final String FORMAT_IDSV_METRICS = "%1$s/%2$s" + COMMON_INITIAL_SUFFIX + ".metrics.idsv.txt";
	private static final String FORMAT_MAPQ_METRICS = "%1$s/%2$s" + COMMON_INITIAL_SUFFIX + ".metrics.mapq.txt";
	private static final String FORMAT_CIGAR_METRICS = "%1$s/%2$s" + COMMON_INITIAL_SUFFIX + ".metrics.cigar.txt";
	private static final String FORMAT_COVERAGE_BLACKLIST_BED = "%1$s/%2$s" + COMMON_INITIAL_SUFFIX + ".coverage.blacklist.bed";
	private static final String FORMAT_REALIGN_FASTQ = "%1$s/%2$s" + COMMON_INITIAL_SUFFIX + ".realign.%3$d.fq";
	private static final String FORMAT_REALIGN_FASTQ_CHR = "%1$s/%3$s/%2$s" + COMMON_INITIAL_SUFFIX + ".%3$s.realign.%4$d.fq";
	private static final String FORMAT_REALIGN_SAM = "%1$s/%2$s" + COMMON_INITIAL_SUFFIX + ".realign.%3$d" + SAM_SUFFIX;
	private static final String FORMAT_REALIGN_SAM_PER_CHR = "%1$s/%3$s/%2$s" + COMMON_INITIAL_SUFFIX + ".%3$s.realign.%4$d" + SAM_SUFFIX;
	private static final String FORMAT_ASSEMBLY_RAW = "%1$s/%2$s" + COMMON_INITIAL_SUFFIX + ".breakend" + SAM_SUFFIX;
	private static final String FORMAT_ASSEMBLY_RAW_PER_CHR = "%1$s/%3$s/%2$s" + COMMON_INITIAL_SUFFIX + ".%3$s.breakend" + SAM_SUFFIX;
	private static final String FORMAT_BREAKPOINT_VCF = "%1$s/%2$s" + COMMON_INITIAL_SUFFIX + ".breakpoint" + VCF_SUFFIX;
	private static final String FORMAT_FLAG_FILE = "%1$s/%2$s" + COMMON_INITIAL_SUFFIX + ".%s";
	private static final String FORMAT_SC_REMOTE_SAM = "%1$s/%2$s" + COMMON_INITIAL_SUFFIX + ".scremote" + SAM_SUFFIX;
	private static final String FORMAT_SC_REMOTE_SAM_PER_CHR = "%1$s/%3$s/%2$s" + COMMON_INITIAL_SUFFIX + ".%3$s.scremote" + SAM_SUFFIX;
	private static final String FORMAT_SC_REMOTE_SAM_UNSORTED = "%1$s/%2$s" + COMMON_INITIAL_SUFFIX + ".scremote.unsorted" + SAM_SUFFIX;
	private static final String FORMAT_SC_REMOTE_SAM_UNSORTED_PER_CHR = "%1$s/%3$s/%2$s" + COMMON_INITIAL_SUFFIX + ".%3$s.scremote.unsorted" + SAM_SUFFIX;;
	private static final String FORMAT_REALIGN_REMOTE_SAM_UNSORTED = "%1$s/%2$s" + COMMON_INITIAL_SUFFIX + ".realignremote.unsorted" + SAM_SUFFIX;
	private static final String FORMAT_REALIGN_REMOTE_SAM_PER_CHR_UNSORTED = "%1$s/%3$s/%2$s" + COMMON_INITIAL_SUFFIX + ".%3$s.realignremote.unsorted" + SAM_SUFFIX;
	private static final String FORMAT_REALIGN_REMOTE_SAM = "%1$s/%2$s" + COMMON_INITIAL_SUFFIX + ".realignremote" + SAM_SUFFIX;
	private static final String FORMAT_REALIGN_REMOTE_SAM_PER_CHR = "%1$s/%3$s/%2$s" + COMMON_INITIAL_SUFFIX + ".%3$s.realignremote" + SAM_SUFFIX;
	private static final String FORMAT_ASSEMBLY = "%1$s/%2$s" + COMMON_INITIAL_SUFFIX + ".assembly" + SAM_SUFFIX;
	private static final String FORMAT_ASSEMBLY_PER_CHR = "%1$s/%3$s/%2$s" + COMMON_INITIAL_SUFFIX + ".%3$s.assembly" + SAM_SUFFIX;
	private static final String FORMAT_ASSEMBLY_MATE = "%1$s/%2$s" + COMMON_INITIAL_SUFFIX + ".assemblymate" + SAM_SUFFIX;
	private static final String FORMAT_ASSEMBLY_MATE_PER_CHR = "%1$s/%3$s/%2$s" + COMMON_INITIAL_SUFFIX + ".%3$s.assemblymate" + SAM_SUFFIX;
	
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
	public synchronized File getIntermediateDirectory(File file) {
		File input = getSource(file);
		File workingDir = new File(getWorkingDirectory(file), input.getName() + INTERMEDIATE_DIR_SUFFIX);
		if (!createdDirectories.contains(workingDir)) {
			if (!workingDir.exists()) {
				workingDir.mkdir();
			}
			createdDirectories.add(workingDir);
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
	private synchronized File getFile(String path) {
		File file = new File(path);
		file.getParentFile().mkdirs();
		return file;
	}
	public File getReadPairBam(File input) {
		return getFile(String.format(FORMAT_RP_SAM, getIntermediateDirectory(input), getSource(input).getName()));
	}
	public File getReadPairBamForChr(File input, String chromosome) {
		return getFile(String.format(FORMAT_RP_SAM_PER_CHR, getIntermediateDirectory(input), getSource(input).getName(), chromosome));
	}
	public File getSoftClipBam(File input) {
		return getFile(String.format(FORMAT_SC_SAM, getIntermediateDirectory(input), getSource(input).getName()));
	}
	public File getSoftClipBamForChr(File input, String chromosome) {
		return getFile(String.format(FORMAT_SC_SAM_PER_CHR, getIntermediateDirectory(input), getSource(input).getName(), chromosome));
	}
	public File getMateBam(File input) {
		return getFile(String.format(FORMAT_MATE_SAM, getIntermediateDirectory(input), getSource(input).getName()));
	}
	public File getMateBamForChr(File input, String chromosome) {
		return getFile(String.format(FORMAT_MATE_SAM_PER_CHR, getIntermediateDirectory(input), getSource(input).getName(), chromosome));
	}
	public File getMateBamUnsorted(File input) {
		return getFile(String.format(FORMAT_MATE_SAM_UNSORTED, getIntermediateDirectory(input), getSource(input).getName()));
	}
	public File getMateBamUnsortedForChr(File input, String chromosome) {
		return getFile(String.format(FORMAT_MATE_SAM_UNSORTED_PER_CHR, getIntermediateDirectory(input), getSource(input).getName(), chromosome));
	}
	public File getAssemblyRawBam(File input) {
		return getFile(String.format(FORMAT_ASSEMBLY_RAW, getIntermediateDirectory(input), getSource(input).getName()));
	}
	public File getAssemblyRawBamForChr(File input, String chromosome) {
		return getFile(String.format(FORMAT_ASSEMBLY_RAW_PER_CHR, getIntermediateDirectory(input), getSource(input).getName(), chromosome));
	}
	public File getAssembly(File input) {
		return getFile(String.format(FORMAT_ASSEMBLY, getIntermediateDirectory(input), getSource(input).getName()));
	}
	public File getAssemblyForChr(File input, String chromosome) {
		return getFile(String.format(FORMAT_ASSEMBLY_PER_CHR, getIntermediateDirectory(input), getSource(input).getName(), chromosome));
	}
	public File getAssemblyMate(File input) {
		return getFile(String.format(FORMAT_ASSEMBLY_MATE, getIntermediateDirectory(input), getSource(input).getName()));
	}
	public File getAssemblyMateForChr(File input, String chromosome) {
		return getFile(String.format(FORMAT_ASSEMBLY_MATE_PER_CHR, getIntermediateDirectory(input), getSource(input).getName(), chromosome));
	}
	public File getBreakpointVcf(File input) {
		return getFile(String.format(FORMAT_BREAKPOINT_VCF, getIntermediateDirectory(input), getSource(input).getName()));
	}
	public File getInsertSizeMetrics(File input) {
		return getFile(String.format(FORMAT_INSERT_SIZE_METRICS, getIntermediateDirectory(input), getSource(input).getName()));
	}
	public File getIdsvMetrics(File input) {
		return getFile(String.format(FORMAT_IDSV_METRICS, getIntermediateDirectory(input), getSource(input).getName()));
	}
	public File getMapqMetrics(File input) {
		return getFile(String.format(FORMAT_MAPQ_METRICS, getIntermediateDirectory(input), getSource(input).getName()));
	}
	public File getCigarMetrics(File input) {
		return getFile(String.format(FORMAT_CIGAR_METRICS, getIntermediateDirectory(input), getSource(input).getName()));
	}
	public File getRealignmentBam(File input, int iteration) {
		return getFile(String.format(FORMAT_REALIGN_SAM, getIntermediateDirectory(input), getSource(input).getName(), iteration));
	}
	public File getRealignmentBamForChr(File input, String chromosome, int iteration) {
		return getFile(String.format(FORMAT_REALIGN_SAM_PER_CHR, getIntermediateDirectory(input), getSource(input).getName(), chromosome, iteration));
	}
	public File getRealignmentFastq(File input, int iteration) {
		return getFile(String.format(FORMAT_REALIGN_FASTQ, getIntermediateDirectory(input), getSource(input).getName(), iteration));
	}
	public File getRealignmentFastqForChr(File input, String chromosome, int iteration) {
		return getFile(String.format(FORMAT_REALIGN_FASTQ_CHR, getIntermediateDirectory(input), getSource(input).getName(), chromosome, iteration));
	}
	public File getFlagFile(File input, String flagname) {
		return getFile(String.format(FORMAT_FLAG_FILE, getIntermediateDirectory(input), getSource(input).getName(), flagname));
	}
	// Remote sorted
	public File getSoftClipRemoteBam(File input) {
		return getFile(String.format(FORMAT_SC_REMOTE_SAM, getIntermediateDirectory(input), getSource(input).getName()));
	}
	public File getSoftClipRemoteBamForChr(File input, String chromosome) {
		return getFile(String.format(FORMAT_SC_REMOTE_SAM_PER_CHR, getIntermediateDirectory(input), getSource(input).getName(), chromosome));
	}
	public File getSoftClipRemoteUnsortedBam(File input) {
		return getFile(String.format(FORMAT_SC_REMOTE_SAM_UNSORTED, getIntermediateDirectory(input), getSource(input).getName()));
	}
	public File getSoftClipRemoteUnsortedBamForChr(File input, String chromosome) {
		return getFile(String.format(FORMAT_SC_REMOTE_SAM_UNSORTED_PER_CHR, getIntermediateDirectory(input), getSource(input).getName(), chromosome));
	}
	public File getRealignmentRemoteUnsortedBam(File input) {
		return getFile(String.format(FORMAT_REALIGN_REMOTE_SAM_UNSORTED, getIntermediateDirectory(input), getSource(input).getName()));
	}
	public File getRealignmentRemoteUnsortedBamForChr(File input, String chromosome) {
		return getFile(String.format(FORMAT_REALIGN_REMOTE_SAM_PER_CHR_UNSORTED, getIntermediateDirectory(input), getSource(input).getName(), chromosome));
	}
	public File getRealignmentRemoteBam(File input) {
		return getFile(String.format(FORMAT_REALIGN_REMOTE_SAM, getIntermediateDirectory(input), getSource(input).getName()));
	}
	public File getRealignmentRemoteBamForChr(File input, String chromosome) {
		return getFile(String.format(FORMAT_REALIGN_REMOTE_SAM_PER_CHR, getIntermediateDirectory(input), getSource(input).getName(), chromosome));
	}
	public File getCoverageBlacklistBed(File input) {
		return getFile(String.format(FORMAT_COVERAGE_BLACKLIST_BED, getIntermediateDirectory(input), getSource(input).getName()));
	}
}
