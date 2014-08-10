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
	private static final String COMMON_INITIAL_SUFFIX = ".idsv";
	private static final String INTERMEDIATE_DIR_SUFFIX = COMMON_INITIAL_SUFFIX + ".working";
	private static final String FORMAT_SC_SAM = "%s" + COMMON_INITIAL_SUFFIX + ".sc.bam";
	private static final String FORMAT_SC_SAM_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s.sc.bam";
	private static final String FORMAT_RP_SAM = "%s" + COMMON_INITIAL_SUFFIX + ".rp.bam";
	private static final String FORMAT_RP_SAM_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s.rp.bam";
	private static final String FORMAT_MATE_SAM = "%s" + COMMON_INITIAL_SUFFIX + ".rpmate.bam";
	private static final String FORMAT_MATE_SAM_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s.rpmate.bam";
	private static final String FORMAT_MATE_SAM_UNSORTED = "%s" + COMMON_INITIAL_SUFFIX + ".rpmate.unsorted.bam";
	private static final String FORMAT_MATE_SAM_UNSORTED_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s.rpmate.unsorted.bam";
	private static final String FORMAT_INSERT_SIZE_METRICS = "%s" + COMMON_INITIAL_SUFFIX + ".metrics.insertsize.txt";
	private static final String FORMAT_IDSV_METRICS = "%s" + COMMON_INITIAL_SUFFIX + ".metrics.idsv.txt";
	private static final String FORMAT_REALIGN_FASTQ = "%s" + COMMON_INITIAL_SUFFIX + ".realign.fq";
	private static final String FORMAT_REALIGN_FASTQ_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s.realign.fq";
	private static final String FORMAT_REALIGN_SAM = "%s" + COMMON_INITIAL_SUFFIX + ".realign.bam";
	private static final String FORMAT_REALIGN_SAM_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s.realign.bam";
	private static final String FORMAT_ASSEMBLY_VCF_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s.breakend.vcf";
	private static final String FORMAT_ASSEMBLY_VCF = "%s" + COMMON_INITIAL_SUFFIX + ".breakend.vcf";
	private static final String FORMAT_BREAKPOINT_VCF_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s-%s.breakpoint.vcf";
	private static final String FORMAT_BREAKPOINT_VCF = "%s" + COMMON_INITIAL_SUFFIX + ".breakpoint.vcf";
	private static final String FORMAT_FLAG_FILE = "%s" + COMMON_INITIAL_SUFFIX + ".%s";
	private static final String FORMAT_SC_REMOTE_SAM = "%s" + COMMON_INITIAL_SUFFIX + ".scremote.bam";
	private static final String FORMAT_SC_REMOTE_SAM_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s.scremote.bam";
	private static final String FORMAT_SC_REMOTE_SAM_UNSORTED = "%s" + COMMON_INITIAL_SUFFIX + ".scremote.unsorted.bam";
	private static final String FORMAT_SC_REMOTE_SAM_UNSORTED_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s.scremote.unsorted.bam";;
	private static final String FORMAT_REALIGN_REMOTE_SAM_UNSORTED = "%s" + COMMON_INITIAL_SUFFIX + ".realignremote.unsorted.bam";
	private static final String FORMAT_REALIGN_REMOTE_SAM_PER_CHR_UNSORTED = "%s" + COMMON_INITIAL_SUFFIX + ".%s.realignremote.unsorted.bam";
	private static final String FORMAT_REALIGN_REMOTE_SAM = "%s" + COMMON_INITIAL_SUFFIX + ".realignremote.bam";
	private static final String FORMAT_REALIGN_REMOTE_SAM_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s.realignremote.bam";
	private static final String FORMAT_ASSEMBLY_REMOTE_VCF = "%s" + COMMON_INITIAL_SUFFIX + ".assemblyremote.vcf";
	private static final String FORMAT_ASSEMBLY_REMOTE_VCF_PER_CHR = "%s" + COMMON_INITIAL_SUFFIX + ".%s.assemblyremote.vcf";
	private static final String FORMAT_ASSEMBLY_REMOTE_VCF_UNSORTED = "%s" + COMMON_INITIAL_SUFFIX + ".assemblyremote.unsorted.vcf";
	private static final String FORMAT_ASSEMBLY_REMOTE_VCF_PER_CHR_UNSORTED = "%s" + COMMON_INITIAL_SUFFIX + ".%s.assemblyremote.unsorted.vcf";
	
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
	public File getAssemblyVcf(File input) {
		return new File(String.format(FORMAT_ASSEMBLY_VCF, getStem(input)));
	}
	public File getAssemblyVcfForChr(File input, String chromosome) {
		return new File(String.format(FORMAT_ASSEMBLY_VCF_PER_CHR, getStem(input), chromosome));
	}
	public File getAssemblyRemoteVcf(File input) {
		return new File(String.format(FORMAT_ASSEMBLY_REMOTE_VCF, getStem(input)));
	}
	public File getAssemblyRemoteVcfForChr(File input, String chromosome) {
		return new File(String.format(FORMAT_ASSEMBLY_REMOTE_VCF_PER_CHR, getStem(input), chromosome));
	}
	public File getAssemblyRemoteUnsortedVcf(File input) {
		return new File(String.format(FORMAT_ASSEMBLY_REMOTE_VCF_UNSORTED, getStem(input)));
	}
	public File getAssemblyRemoteUnsortedVcfForChr(File input, String chromosome) {
		return new File(String.format(FORMAT_ASSEMBLY_REMOTE_VCF_PER_CHR_UNSORTED, getStem(input), chromosome));
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
