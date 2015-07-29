package au.edu.wehi.idsv;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMSequenceRecord;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public abstract class EvidenceSource {
	public abstract int getMaxConcordantFragmentSize();
	public abstract int getMinConcordantFragmentSize();
	protected final File input;
	private final ProcessingContext processContext;
	public abstract int getRealignmentIterationCount();
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
		for (int i = 0; i < getRealignmentIterationCount(); i++) {
			if (processContext.shouldProcessPerChromosome()) {
				for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
					sb.append(getBowtie2Script(threads, fsc.getRealignmentFastqForChr(input, seq.getSequenceName(), i), fsc.getRealignmentBamForChr(input, seq.getSequenceName(), i)));
				}
			} else {
				sb.append(getBowtie2Script(threads, fsc.getRealignmentFastq(input, i), fsc.getRealignmentBam(input, i)));
			}
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
		List<File> bamList = new ArrayList<File>();
		List<File> fastqList = new ArrayList<File>();
		for (int i = 0; i < getRealignmentIterationCount(); i++) {
			if (processContext.shouldProcessPerChromosome()) {
				for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
					bamList.add(fsc.getRealignmentBamForChr(input, seq.getSequenceName(), i));
					fastqList.add(fsc.getRealignmentFastqForChr(input, seq.getSequenceName(), i));
				}
			} else {
				bamList.add(fsc.getRealignmentBam(input, i));
				fastqList.add(fsc.getRealignmentFastq(input, i));
			}
		}
		return isRealignmentComplete(bamList, fastqList);
	}
	private boolean isRealignmentComplete(List<File> bamList, List<File> fastqList) {
		assert(bamList.size() == fastqList.size());
		boolean complete = IntermediateFileUtil.checkIntermediate(bamList, fastqList);
		if (!complete) {
			// mock up empty files for the chromosomes we have no reads to realign
			for (int i = 0; i < bamList.size(); i++) {
				File bam = bamList.get(i);
				File fastq = fastqList.get(i);
				if (fastq.exists() && fastq.length() == 0) {
					if (!bam.exists() || bam.lastModified() < fastq.lastModified()) {
						writeEmptyBAM(bam);
					}
				}
			}
			complete = IntermediateFileUtil.checkIntermediate(bamList, fastqList);
		}
		return complete; 
	}
	private void writeEmptyBAM(File sam) {
		SAMFileWriter w = processContext.getSamFileWriterFactory(false).makeSAMOrBAMWriter(processContext.getBasicSamHeader(), false, sam);
		w.close();
	}
}