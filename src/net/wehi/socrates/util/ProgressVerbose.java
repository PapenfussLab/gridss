/**
 * 
 */
package net.wehi.socrates.util;


/**
 * @author hsu
 *
 * Created on Jan 23, 2013
 */
public class ProgressVerbose {
	private long interval = 0;
	private long progress = 0;
	private boolean verbose = true;
	
	public ProgressVerbose(String initMessage, long verboseInterval, boolean verbose) {
		if (verbose) System.out.println( initMessage );
		interval = verboseInterval;
		this.verbose = verbose;
	}
	
	public void stepProgress(String suffixMessage) {
		progress++;
		if (this.verbose && progress % interval == 0) {
			System.err.println( progress + " " + suffixMessage);
		}
	}
	
	public void end(String suffixMessage) {
		if (this.verbose) System.err.println( progress + " " + suffixMessage);
	}
	
	public void stepProgress(String prefixMessage, String suffixMessage) {
		progress++;
		if (this.verbose && progress % interval == 0) {
			System.err.println( prefixMessage + " " + progress + " " + suffixMessage);
		}
	}
	
	public void end(String prefixMessage, String suffixMessage) {
		if (this.verbose) System.err.println( prefixMessage + " " + progress + " " + suffixMessage);
	}

}
