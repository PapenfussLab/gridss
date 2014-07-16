package au.edu.wehi.idsv;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Log;

import java.io.File;

/**
 * Utility helper class for intermediate files
 * @author cameron.d
 *
 */
public abstract class IntermediateFileUtil {
	private static final Log log = Log.getInstance(IntermediateFileUtil.class);
	/**
	 * Checks that the given intermediate file is valid
	 * @param file file to check
	 * @param source source file intermediate has been generated from
	 * @return true if intermediate file appears to be valid
	 */
	public static boolean checkIntermediate(File file, File source) {
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
	public static  boolean checkIntermediate(File file) {
		return checkIntermediate(file, null);
	}
}
