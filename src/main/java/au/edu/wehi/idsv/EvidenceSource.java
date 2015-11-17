package au.edu.wehi.idsv;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.StringUtils;

import au.edu.wehi.idsv.alignment.ExternalProcessFastqAligner;
import au.edu.wehi.idsv.alignment.FastqAligner;

import com.google.common.io.Files;

public abstract class EvidenceSource {
	private static final Log log = Log.getInstance(EvidenceSource.class);
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
	public String getRealignmentScript() {
		if (isRealignmentComplete(true)) return "";
		return getBowtie2Script();
	}
	private String getBowtie2Script() {
		StringBuilder sb = new StringBuilder();
		FileSystemContext fsc = processContext.getFileSystemContext();
		for (int i = 0; i < getRealignmentIterationCount(); i++) {
			if (processContext.shouldProcessPerChromosome()) {
				for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
					sb.append(getBowtie2Script(fsc.getRealignmentFastqForChr(input, seq.getSequenceName(), i), fsc.getRealignmentBamForChr(input, seq.getSequenceName(), i)));
				}
			} else {
				sb.append(getBowtie2Script(fsc.getRealignmentFastq(input, i), fsc.getRealignmentBam(input, i)));
			}
		}
		return sb.toString();
	}
	private String getBowtie2Script(File realignFastq, File realignBam) {
		if (realignFastq.exists() && !realignBam.exists()) {
			return String.format("bowtie2 --threads %d --local --mm --reorder -x \"%s\" -U \"%s\" | samtools view -Sb -o \"%s\" - \n",
					processContext.getWorkerThreadCount(),
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
	public boolean isRealignmentComplete(boolean performRealignment) {
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
					if (performRealignment) {
						if (i > 0) {
							generateFastqFromBam(fsc.getRealignmentBamForChr(input, seq.getSequenceName(), i - 1), fastq);
						}
						align(fastq, bam);	
					}
				}
			} else {
				File bam = fsc.getRealignmentBam(input, i);
				File fastq = fsc.getRealignmentFastq(input, i);
				bamList.add(bam);
				fastqList.add(fastq);
				if (performRealignment) {
					if (i > 0) {
						generateFastqFromBam(fsc.getRealignmentBam(input, i - 1), fastq);
					}
					align(fastq, bam);	
				}
			}
		}
		return IntermediateFileUtil.checkIntermediate(bamList, fastqList, getContext().getConfig().ignoreFileTimestamps);
	}
	/**
	 * Generates a compound breakpoint fastq from the previous partial alignments 
	 * @param bam previous breakpoint alignments
	 * @param fastq soft clips alignments requiring realignment for compound breakpoint identification
	 */
	private void generateFastqFromBam(File bam, File fastq) {
		if (!bam.exists()) return;
		if (!IntermediateFileUtil.checkIntermediate(fastq, bam, getContext().getConfig().ignoreFileTimestamps)) {
			// need to create the fastq
			try (SamReader reader = processContext.getSamReader(bam)) {
				try (CloseableIterator<SAMRecord> it = processContext.getSamReaderIterator(reader)) {
					File tmp = FileSystemContext.getWorkingFileFor(fastq);
					try (FastqBreakpointWriter writer = new FastqBreakpointWriter(processContext.getFastqWriterFactory().newWriter(tmp))) {
						while (it.hasNext()) {
							SAMRecord read = it.next();
							writer.writeRealignment(read, processContext.getRealignmentParameters().minLength);
						}
					}
					Files.move(tmp, fastq);
					log.debug("Generated ", fastq);
				}
			} catch (IOException e) {
				log.error(e, "Error generating " + fastq.toString());
				throw new RuntimeException(e);
			}
		}
	}
	private void align(File fastq, File bam) {
		if (IntermediateFileUtil.checkIntermediate(bam, fastq, getContext().getConfig().ignoreFileTimestamps)) {
			// bam exists - no need to realign
			return;
		}
		if (fastq.exists() && fastq.length() == 0) {
			// no need to make call to external aligner when we have 0 record to align
			writeEmptyBAM(bam);
			return;
		}
		List<String> cmdline = processContext.getConfig().getRealignment().commandline;
		if (cmdline != null && cmdline.size() > 0 && !StringUtils.isEmpty(cmdline.get(0))) {
			File tmp = FileSystemContext.getWorkingFileFor(bam, "aligning.");
			FastqAligner aligner = new ExternalProcessFastqAligner(processContext.getSamReaderFactory(), processContext.getSamFileWriterFactory(false), cmdline);
			try {
				aligner.align(fastq, tmp, processContext.getReferenceFile(), processContext.getWorkerThreadCount());
				Files.move(tmp, bam);
			} catch (IOException e) {
				log.error(e);
				throw new RuntimeException(e);
			}
		}
	}
	private void writeEmptyBAM(File sam) {
		try (SAMFileWriter w = processContext.getSamFileWriterFactory(false).makeSAMOrBAMWriter(processContext.getBasicSamHeader(), false, sam)) {
		}
	}
}