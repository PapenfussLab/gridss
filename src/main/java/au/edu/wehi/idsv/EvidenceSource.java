package au.edu.wehi.idsv;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMSequenceRecord;

import java.io.File;

public abstract class EvidenceSource {
	public abstract int getMaxConcordantFragmentSize();
	protected final File input;
	private final ProcessingContext processContext;
	public ProcessingContext getContext() {
		return processContext;
	}
	/**
	 * Gets the file that the intermediate directory location and stucture is based on.
	 * @return anchor file
	 */
	public File getFileIntermediateDirectoryBasedOn() { return input; }
	/**
	 * New evidence source
	 * @param input base file for which intermediates files are relative to
	 */
	public EvidenceSource(ProcessingContext processContext, File input) {
		this.processContext = processContext;
		this.input = input;
	}
	public String getRealignmentScript(int threads) {
		return getBowtie2Script(threads);
	}
	private String getBowtie2Script(int threads) {
		StringBuilder sb = new StringBuilder();
		FileSystemContext fsc = processContext.getFileSystemContext();
		if (processContext.shouldProcessPerChromosome()) {
			for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
				sb.append(getBowtie2Script(threads, fsc.getRealignmentFastqForChr(input, seq.getSequenceName()), fsc.getRealignmentBamForChr(input, seq.getSequenceName())));
			}
		} else {
			sb.append(getBowtie2Script(threads, fsc.getRealignmentFastq(input), fsc.getRealignmentBam(input)));
		}
		return sb.toString();
	}
	private String getBowtie2Script(int threads, File realignFastq, File realignBam) {
		if (realignFastq.exists() && realignFastq.length() > 0 && (!realignBam.exists() || realignFastq.lastModified() > realignBam.lastModified())) {
			return String.format("bowtie2 --threads %d --local --mm --reorder -x \"%s\" -U \"%s\" | samtools view -Sb -o \"%s\" - \n",
					threads,
					processContext.getReferenceFile(),
					realignFastq,
					realignBam);
		}
		return "";
	}
	/**
	 * Checks if realignment is complete for the given source file
	 * @param processContext
	 * @param source source
	 * @return
	 */
	public boolean isRealignmentComplete() {
		FileSystemContext fsc = processContext.getFileSystemContext();
		if (processContext.shouldProcessPerChromosome()) {
			for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
				if (!isRealignmentComplete(fsc.getRealignmentBamForChr(input, seq.getSequenceName()), fsc.getRealignmentFastqForChr(input, seq.getSequenceName()))) {
					return false;
				}
			}
		} else {
			if (!isRealignmentComplete(fsc.getRealignmentBam(input), fsc.getRealignmentFastq(input))) {
				return false;
			}
		}
		return true;
	}
	private boolean isRealignmentComplete(File bam, File fastq) {
		if (fastq.exists() && fastq.length() == 0) {
			if (!bam.exists() || bam.lastModified() < fastq.lastModified()) {
				writeEmptyBAM(bam);
			}
		}
		return IntermediateFileUtil.checkIntermediate(bam, fastq);
	}
	private void writeEmptyBAM(File sam) {
		SAMFileWriter w = processContext.getSamFileWriterFactory(false).makeSAMOrBAMWriter(processContext.getBasicSamHeader(), false, sam);
		w.close();
	}
}