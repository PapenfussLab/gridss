package au.edu.wehi.idsv;

import htsjdk.samtools.SAMSequenceRecord;

import java.io.File;

public abstract class EvidenceSource {
	public abstract int getMaxConcordantFragmentSize();
	protected final File input;
	protected final ProcessingContext processContext;
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
		return String.format("bowtie2 --threads %d --local --mm --reorder -x \"%s\" -U \"%s\" | samtools view -Sb -o \"%s\" - \n",
				threads,
				processContext.getReferenceFile(),
				realignFastq,
				realignBam);
	}
	/**
	 * Checks if realignment is complete for the given source file
	 * @param processContext
	 * @param source source
	 * @return
	 */
	public boolean isRealignmentComplete() {
		boolean done = true;
		FileSystemContext fsc = processContext.getFileSystemContext();
		if (processContext.shouldProcessPerChromosome()) {
			for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
				done &= IntermediateFileUtil.checkIntermediate(fsc.getRealignmentBamForChr(input, seq.getSequenceName()), fsc.getRealignmentFastqForChr(input, seq.getSequenceName()));
				if (!done) return false;
			}
		} else {
			done &= IntermediateFileUtil.checkIntermediate(fsc.getRealignmentBam(input), fsc.getRealignmentFastq(input));
			if (!done) return false;
		}
		return done;
	}
}