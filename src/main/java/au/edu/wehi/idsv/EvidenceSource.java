package au.edu.wehi.idsv;

import java.io.File;

public abstract class EvidenceSource {
	public abstract int getMaxConcordantFragmentSize();
	public abstract int getMinConcordantFragmentSize();
	protected final File file;
	private final ProcessingContext processContext;
	public ProcessingContext getContext() {
		return processContext;
	}
	/**
	 * Gets the file that the intermediate directory location and stucture is based on.
	 * @return anchor file
	 */
	public File getFile() { return file; }
	/**
	 * New evidence source
	 * @param input base file for which intermediates files are relative to
	 */
	public EvidenceSource(ProcessingContext processContext, File file) {
		this.processContext = processContext;
		this.file = file;
	}
}