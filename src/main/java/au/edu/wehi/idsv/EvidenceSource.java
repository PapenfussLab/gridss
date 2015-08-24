package au.edu.wehi.idsv;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;

import java.io.File;
import java.io.IOException;
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
		if (isRealignmentComplete()) return "";
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
		if (realignFastq.exists() && !realignBam.exists()) {
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
					File bam = fsc.getRealignmentBamForChr(input, seq.getSequenceName(), i);
					File fastq = fsc.getRealignmentFastqForChr(input, seq.getSequenceName(), i);
					bamList.add(bam);
					fastqList.add(fastq);
					alignIfEmpty(fastq, bam);
				}
			} else {
				File bam = fsc.getRealignmentBam(input, i);
				File fastq = fsc.getRealignmentFastq(input, i);
				bamList.add(bam);
				fastqList.add(fastq);
				alignIfEmpty(fastq, bam);
			}
		}
		return IntermediateFileUtil.checkIntermediate(bamList, fastqList);
	}
	public void iterateRealignment() throws IOException {
		FileSystemContext fsc = processContext.getFileSystemContext();
		for (int i = 1; i < getRealignmentIterationCount(); i++) {
			if (processContext.shouldProcessPerChromosome()) {
				for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
					File lastbam = fsc.getRealignmentBamForChr(input, seq.getSequenceName(), i - 1);
					File lastfastq = fsc.getRealignmentFastqForChr(input, seq.getSequenceName(), i - 1);
					File bam = fsc.getRealignmentBamForChr(input, seq.getSequenceName(), i);
					File fastq = fsc.getRealignmentFastqForChr(input, seq.getSequenceName(), i);
					iterateRealignment(lastbam, lastfastq, bam, fastq);
				}
			} else {
				File lastbam = fsc.getRealignmentBam(input, i - 1);
				File lastfastq = fsc.getRealignmentFastq(input, i - 1);
				File bam = fsc.getRealignmentBam(input, i);
				File fastq = fsc.getRealignmentFastq(input, i);
				iterateRealignment(lastbam, lastfastq, bam, fastq);
			}
		}
	}
	private void iterateRealignment(File lastbam, File lastfastq, File bam, File fastq) throws IOException {
		if (bam.exists()) return;
		if (fastq.exists()) return;
		if (!lastbam.exists()) return;
		// need to create the fastq for the next iteration of realignment
		SamReader reader = processContext.getSamReader(lastbam);
		CloseableIterator<SAMRecord> it = processContext.getSamReaderIterator(reader);
		FastqBreakpointWriter writer = new FastqBreakpointWriter(processContext.getFastqWriterFactory().newWriter(fastq));
		while (it.hasNext()) {
			SAMRecord read = it.next();
			writer.writeRealignment(read, processContext.getRealignmentParameters().minLength);
		}		
		it.close();
		reader.close();
		writer.close();
		alignIfEmpty(fastq, bam);
	}
	private void alignIfEmpty(File fastq, File bam) {
		if (fastq.exists() && fastq.length() == 0 && !bam.exists()) {
			// no need to make call to external aligner when we have 0 record to align
			writeEmptyBAM(bam);
		}
	}
	private void writeEmptyBAM(File sam) {
		SAMFileWriter w = processContext.getSamFileWriterFactory(false).makeSAMOrBAMWriter(processContext.getBasicSamHeader(), false, sam);
		w.close();
	}
}