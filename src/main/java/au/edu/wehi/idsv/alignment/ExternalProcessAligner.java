package au.edu.wehi.idsv.alignment;

import java.io.File;
import java.util.concurrent.TimeUnit;

import htsjdk.samtools.util.Log;

class ExternalProcessHelper {
	/**
	 * Number of seconds to wait for the external process to shut down
	 */
	private static final int SHUTDOWN_GRACE_PERIOD = 60;
	private static final Log log = Log.getInstance(ExternalProcessHelper.class);	
	public static void shutdown(Process process, String commandline, String helpfulErrorMessage) {
		try {
			process.waitFor(SHUTDOWN_GRACE_PERIOD, TimeUnit.SECONDS);
		} catch (InterruptedException e) {
			// Restore the interrupted status
            Thread.currentThread().interrupt();
		}
		if (process.isAlive()) {
			log.error(String.format("External process still alive: \"%s\"", commandline));
			try {
				process.destroyForcibly().waitFor(SHUTDOWN_GRACE_PERIOD, TimeUnit.SECONDS);
			} catch (InterruptedException e) {
				Thread.currentThread().interrupt();
			}
			if (process.isAlive()) {
				log.error(String.format("External process still alive after being forcibly terminated: \"%s\"", commandline));
			}
		}
		if (process.exitValue() != 0) {
			String msg = String.format(
					"Subprocess terminated with with exit status %d. "
					+ "Failed executing \"%s\". %s",
					process.exitValue(),
					commandline,
					helpfulErrorMessage);
			log.error(msg);
			throw new RuntimeException(msg);
		}
	}
	public static void shutdownAligner(Process process, String commandline, File reference) {
		shutdown(process, commandline, String.format(
				"Can you run the alignment command from the command line? "
				+ "Is the aligner on PATH? "
				+ "Did you build an index with prefix %s?",
				reference));
	}
}
