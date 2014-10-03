package au.edu.wehi.idsv.pipeline;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;

import java.io.Closeable;
import java.io.File;
import java.util.EnumSet;
import java.util.List;

import au.edu.wehi.idsv.ProcessStep;
import au.edu.wehi.idsv.ProcessingContext;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

public abstract class DataTransformStep implements Closeable {
	protected abstract Log getLog();
	protected final List<Object> toClose = Lists.newArrayList();
	protected final ProcessingContext processContext;
	public DataTransformStep(ProcessingContext processContext) {
		this.processContext = processContext;
	}
	/**
	 * Outputs from performing step
	 * @return output files
	 */
	public abstract List<File> getOutput();
	/**
	 * Temporary files used during processing
	 * @return temporary files
	 */
	public List<File> getTemporary() { return ImmutableList.of(); }
	/**
	 * Input files used for processing
	 * @return input files
	 */
	public abstract List<File> getInputs();
	/**
	 * Determines whether this transformation has been completed
	 * @return true if complete, false otherwise
	 */
	public boolean isComplete() {
		return allExist(getOutput());
	}
	/**
	 * Determines if this transformation is able to run
	 * @return true if able to run, false otherwise
	 */
	public boolean canProcess() {
		return allExist(getInputs());
	}
	private boolean allExist(Iterable<File> files) {
		boolean exists = true; 
		for (File f : files) {
			exists &= f.exists();
		}
		return exists;
	}
	private void tryDelete(File f) {
		try {
			if (f.exists()) {
				if (!f.delete()) {
					getLog().error("Unable to delete intermediate file ", f,  " during rollback.");
				}
			}
		} catch (Exception e) {
			getLog().error(e, "Unable to delete intermediate file ", f,  " during rollback.");
		}
	}
	/**
	 * Deletes all temporary files
	 */
	protected void deleteTemp() {
		for (File f : getTemporary()) {
			tryDelete(f);
		}
	}
	/**
	 * Deletes all output files
	 */
	protected void deleteOutput() {
		for (File f : getOutput()) {
			tryDelete(f);
		}
	}
	public void close() {
		for (Object obj : toClose) {
			CloserUtil.close(obj);
		}
		toClose.clear();
	}
	public abstract EnumSet<ProcessStep> process(EnumSet<ProcessStep> steps);
}
