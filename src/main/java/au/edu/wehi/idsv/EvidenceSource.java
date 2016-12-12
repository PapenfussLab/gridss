package au.edu.wehi.idsv;

import java.io.File;

import htsjdk.samtools.SAMFileHeader.SortOrder;

public abstract class EvidenceSource {
	public abstract int getMaxConcordantFragmentSize();
	public abstract int getMinConcordantFragmentSize();
	public abstract int getMaxReadLength();
	public abstract int getMaxReadMappedLength();
	private final ProcessingContext processContext;
	private final File file;
	private final File nameSorted;
	public ProcessingContext getContext() {
		return processContext;
	}
	/**
	 * Gets the file that the intermediate directory location and stucture is based on.
	 * @return anchor file
	 */
	public File getFile() { return file; }
	/**
	 * Gets the file that the intermediate directory location and stucture is based on.
	 * @return anchor file
	 */
	public File getFile(SortOrder sortOrder) {
		switch (sortOrder) {
			case coordinate:
			case unsorted:
				return getFile();
			case queryname:
				return nameSorted;
			default:
				throw new IllegalArgumentException(String.format("Unhandled sort order %s", sortOrder));
		}
	}
	/**
	 * New evidence source
	 * @param nameSorted TODO
	 * @param input base file for which intermediates files are relative to
	 */
	public EvidenceSource(ProcessingContext processContext, File file, File nameSorted) {
		this.processContext = processContext;
		this.file = file;
		this.nameSorted = nameSorted;
	}
}