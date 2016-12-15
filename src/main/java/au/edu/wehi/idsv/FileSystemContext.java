package au.edu.wehi.idsv;

import java.io.File;
import java.util.HashSet;
import java.util.Set;

import gridss.analysis.CollectCigarMetrics;
import gridss.analysis.CollectIdsvMetrics;
import gridss.analysis.CollectMapqMetrics;
import gridss.analysis.CollectTagMetrics;

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
		return getWorkingFileFor(file, "gridss.tmp.");
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
	private static final String COMMON_INITIAL_SUFFIX = ".gridss";
	private static final String INTERMEDIATE_DIR_SUFFIX = COMMON_INITIAL_SUFFIX + ".working";
	private static final String FORMAT_SV_SAM = "%1$s/%2$s.sv.bam";
	private static final String FORMAT_METRICS_PREFIX = "%1$s/%2$s";
	private static final String FORMAT_INSERT_SIZE_METRICS = FORMAT_METRICS_PREFIX + ".insert_size_metrics";
	private static final String FORMAT_IDSV_METRICS = FORMAT_METRICS_PREFIX + CollectIdsvMetrics.METRICS_SUFFIX;
	private static final String FORMAT_MAPQ_METRICS = FORMAT_METRICS_PREFIX + CollectMapqMetrics.METRICS_SUFFIX;
	private static final String FORMAT_CIGAR_METRICS = FORMAT_METRICS_PREFIX + CollectCigarMetrics.METRICS_SUFFIX;
	private static final String FORMAT_TAG_METRICS = FORMAT_METRICS_PREFIX + CollectTagMetrics.METRICS_SUFFIX;
	private static final String FORMAT_COVERAGE_BLACKLIST_BED = "%1$s/%2$s.coverage.blacklist.bed";
	private static final String FORMAT_REALIGN_FASTQ = "%1$s/%2$s.realign.%3$d.fq";
	private static final String FORMAT_REALIGN_SAM = "%1$s/%2$s.realign.%3$d" + SAM_SUFFIX;
	private static final String FORMAT_BREAKPOINT_VCF = "%1$s/%2$s.breakpoint" + VCF_SUFFIX;
	private static final String FORMAT_ASSEMBLY_CHUNK_SAM = "%1$s/%2$s.assembly.chunk%3$d" + SAM_SUFFIX;
	
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
	public File getSVBam(File input) {
		return getFile(String.format(FORMAT_SV_SAM, getIntermediateDirectory(input), getSource(input).getName()));
	}
	public File getBreakpointVcf(File input) {
		return getFile(String.format(FORMAT_BREAKPOINT_VCF, getIntermediateDirectory(input), getSource(input).getName()));
	}
	public File getMetricsPrefix(File input) {
		return getFile(String.format(FORMAT_METRICS_PREFIX, getIntermediateDirectory(input), getSource(input).getName()));
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
	public File getTagMetrics(File input) {
		return getFile(String.format(FORMAT_TAG_METRICS, getIntermediateDirectory(input), getSource(input).getName()));
	}
	public File getRealignmentBam(File input, int iteration) {
		return getFile(String.format(FORMAT_REALIGN_SAM, getIntermediateDirectory(input), getSource(input).getName(), iteration));
	}
	public File getRealignmentFastq(File input, int iteration) {
		return getFile(String.format(FORMAT_REALIGN_FASTQ, getIntermediateDirectory(input), getSource(input).getName(), iteration));
	}
	public File getCoverageBlacklistBed(File input) {
		return getFile(String.format(FORMAT_COVERAGE_BLACKLIST_BED, getIntermediateDirectory(input), getSource(input).getName()));
	}
	public File getAssemblyChunkBam(File input, int chunk) {
		return getFile(String.format(FORMAT_ASSEMBLY_CHUNK_SAM, getIntermediateDirectory(input), getSource(input).getName(), chunk));
	}
	
}
