package au.edu.wehi.idsv;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.util.NoSuchElementException;

import au.edu.wehi.idsv.metrics.RelevantMetrics;

import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

public abstract class EvidenceSource implements Iterable<DirectedEvidence> {
	protected abstract boolean isProcessingComplete();
	protected abstract void process();
	protected abstract CloseableIterator<DirectedEvidence> perChrIterator(String chr);
	protected abstract CloseableIterator<DirectedEvidence> singleFileIterator();
	public abstract RelevantMetrics getMetrics();
	private static final Log log = Log.getInstance(EvidenceSource.class);
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
	/**
	 * Checks that the given intermediate file is valid
	 * @param file file to check
	 * @param source source file
	 * @return true if intermediate file appears to be valid
	 */
	protected boolean checkIntermediate(File file, File source) {
		if (!file.exists()) {
			log.debug("Missing intermediate ", file);
			return false;
		}
		if (source != null && source.exists() && file.lastModified() < source.lastModified()) {
			log.info(source, " has a more recent timestamp than ", file, ". Considering ", file, " out of date.");
			return false;
		}
		return true;
	}
	/**
	 * Checks that the given intermediate file is valid
	 * @param file file to check
	 * @return true if intermediate file appears to be valid
	 */
	protected boolean checkIntermediate(File file) {
		return checkIntermediate(file, null);
	}
	public boolean isRealignmentComplete() {
		boolean done = true;
		FileSystemContext fsc = processContext.getFileSystemContext();
		if (processContext.shouldProcessPerChromosome()) {
			for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
				done &= checkIntermediate(fsc.getRealignmentBamForChr(input, seq.getSequenceName()), fsc.getSVBamForChr(input, seq.getSequenceName()));
			}
		} else {
			done &= checkIntermediate(fsc.getRealignmentBam(input), fsc.getSVBam(input));
		}
		return done;
	}
	public String getRealignmentScript() {
		return getBowtie2Script();
	}
	private String getBowtie2Script() {
		StringBuilder sb = new StringBuilder();
		FileSystemContext fsc = processContext.getFileSystemContext();
		if (processContext.shouldProcessPerChromosome()) {
			for (SAMSequenceRecord seq : processContext.getReference().getSequenceDictionary().getSequences()) {
				sb.append(getBowtie2Script(fsc.getRealignmentFastqForChr(input, seq.getSequenceName()), fsc.getRealignmentBamForChr(input, seq.getSequenceName())));
			}
		} else {
			sb.append(getBowtie2Script(fsc.getRealignmentFastq(input), fsc.getRealignmentBam(input)));
		}
		return sb.toString();
	}
	private String getBowtie2Script(File realignFastq, File realignBam) {
		return String.format("bowtie2 --local --mm --reorder -x \"%s\" -U \"%s\" | samtools view -b -o \"%s\" - \n",
				processContext.getReferenceFile(),
				realignFastq,
				realignBam);
	}
	public CloseableIterator<DirectedEvidence> iterator() {
		if (!isProcessingComplete()) {
			process();
		}
		if (!isProcessingComplete()) {
			throw new IllegalStateException(String.format("Missing intermediate files for ", input));
		} else if (!isRealignmentComplete()) {
			log.debug("Missing realignment bam. Traversing breakends only.");
		}
		if (processContext.shouldProcessPerChromosome()) {
			// Lazily iterator over each input
			return new PerChromosomeAggregateIterator();
		} else {
			return singleFileIterator();
		}
	}
	public CloseableIterator<DirectedEvidence> iterator(String chr) {
		if (!isProcessingComplete()) {
			process();
		}
		if (!isProcessingComplete()) {
			throw new IllegalStateException(String.format("Missing intermediate files for ", input));
		} else if (!isRealignmentComplete()) {
			log.debug("Missing realignment bam. Traversing breakends only.");
		}
		if (processContext.shouldProcessPerChromosome()) {
			return perChrIterator(chr);
		} else {
			return new CloseableChromosomeFilterIterator(singleFileIterator(), processContext.getDictionary().getSequence(chr).getSequenceIndex());
		}
	}
	/**
	 * Lazy iterator that iterates over all per chromosome evidence
	 * only opening new file when required and closing as soon as possible 
	 * @author Daniel Cameron
	 *
	 */
	private class PerChromosomeAggregateIterator implements CloseableIterator<DirectedEvidence> {
		private int currentReferenceIndex = -1;
		private CloseableIterator<DirectedEvidence> currentSource = null;
		private boolean hasMoreContigs() {
			return currentReferenceIndex < processContext.getReference().getSequenceDictionary().size() - 1;
		}
		private boolean currentHasNext() {
			if (currentSource != null) {
				if (currentSource.hasNext()) {
					return true;
				} else {
					// close as soon as we know there's no more records left
					closeCurrent();
				}
			}
			return false;
		}
		private boolean advanceContig() {
			assert(hasMoreContigs());
			closeCurrent();
			currentReferenceIndex++;
			currentSource = perChrIterator(processContext.getReference().getSequenceDictionary().getSequence(currentReferenceIndex).getSequenceName());
			return true;
		}
		private void closeCurrent() {
			if (currentSource != null) {
				currentSource.close();
			}
			currentSource = null;
		}
		@Override
		public boolean hasNext() {
			while (!currentHasNext() && hasMoreContigs()) {
				advanceContig();
			}
			return currentHasNext();
		}
		@Override
		public DirectedEvidence next() {
			if (hasNext()) return currentSource.next();
			throw new NoSuchElementException();
		}
		@Override
		public void close() {
			closeCurrent();
			currentReferenceIndex = processContext.getReference().getSequenceDictionary().size();
		}
	}
	/**
	 * Filters the iterable to only results on the given chromsome
	 * closing the underlying resources as soon as possible
	 * 
	 * @author Daniel Cameron
	 *
	 */
	private static class CloseableChromosomeFilterIterator implements CloseableIterator<DirectedEvidence> {
		private CloseableIterator<DirectedEvidence> source;
		private PeekingIterator<DirectedEvidence> it;
		private int referenceIndex;
		private boolean closed = false;
		public CloseableChromosomeFilterIterator(CloseableIterator<DirectedEvidence> it, int referenceIndex) {
			assert(referenceIndex >= 0);
			this.source = it;
			this.it = Iterators.peekingIterator(it);
			this.referenceIndex = referenceIndex;
		}
		@Override
		public boolean hasNext() {
			if (closed) return false;
			while (it.hasNext() && it.peek().getBreakendSummary().referenceIndex < referenceIndex) {
				// fast-forward advance to our chr
				it.next();
			}
			if (it.hasNext() && it.peek().getBreakendSummary().referenceIndex > referenceIndex) {
				// if we've nothing left on our chr, close immediately
				close();
				return false;
			}
			if (!it.hasNext()) {
				close();
				return false;
			}
			return !closed && it.hasNext();
		}
		@Override
		public DirectedEvidence next() {
			if (hasNext()) return it.next();
			throw new NoSuchElementException();
		}
		@Override
		public void close() {
			if (!closed) {
				source.close();
				source = null;
				it = null;
			}
			closed = true;
		}
	}
}