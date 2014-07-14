package au.edu.wehi.idsv;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.util.Iterator;

public abstract class EvidenceSource implements Iterable<DirectedEvidence> {
	public abstract Iterator<DirectedEvidence> iterator(String chr);
	private static final Log log = Log.getInstance(EvidenceSource.class);
	protected final File input;
	protected final ProcessingContext processContext;
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
				done |= checkIntermediate(fsc.getRealignmentBamForChr(input, seq.getSequenceName()), fsc.getSVBamForChr(input, seq.getSequenceName()));
			}
		} else {
			done |= checkIntermediate(fsc.getRealignmentBam(input), fsc.getSVBam(input));
		}
		return done;
	}
	public String getRealignmentScript() {
		
	}
	public String getBowtie2Script() {
		
	}
}
