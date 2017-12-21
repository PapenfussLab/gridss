package au.edu.wehi.idsv;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import htsjdk.samtools.util.Log;

/**
 * Utility helper class for intermediate files
 * @author Daniel Cameron
 *
 */
public abstract class IntermediateFileUtil {
	private static final Log log = Log.getInstance(IntermediateFileUtil.class);
	/**
	 * Timestamps this close together are consider written at the same time.
	 */
	private static final int TIMESTAMP_MS_ALLOWABLE_ERROR = 0;
	private static SimpleDateFormat getLogDateFormat() { return new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.SSS"); }
	/**
	 * Checks that the given intermediate file is valid
	 * @param file file to check
	 * @param source source file intermediate has been generated from
	 * @return true if intermediate file appears to be valid
	 */
	public static boolean checkIntermediate(File file, File source) {
		if (!file.exists()) {
			//log.debug("Missing intermediate ", file);
			return false;
		}
		if (!gridss.Defaults.IGNORE_FILE_TIMESTAMPS && source != null && source.exists() && file.lastModified() < source.lastModified() - TIMESTAMP_MS_ALLOWABLE_ERROR) {
			log.info(source, " has a more recent timestamp (",
					getLogDateFormat().format(new Date(source.lastModified())),") than ",
					file, " (",
					getLogDateFormat().format(new Date(file.lastModified())),
					"). Considering the latter out of date.");
			return false;
		}
		return true;
	}
	/**
	 * Checks that the given intermediate file is valid
	 * @param file file to check
	 * @return true if intermediate file appears to be valid
	 */
	public static boolean checkIntermediate(File file) {
		return checkIntermediate(file, null);
	}
	public static boolean checkIntermediate(List<File> fileList, List<File> sourceList) {
		assert(fileList.size() == sourceList.size());
		Map<String, File[]> dirList = new HashMap<String, File[]>();
		for (int i = 0; i < fileList.size(); i++) {
			File file = fileList.get(i);
			File source = sourceList.get(i);
			if (!checkIntermediateCached(dirList, file, source)) {
				return false;
			}
		}
		return true;
	}
	private static boolean checkIntermediateCached(Map<String, File[]> dirList, File file, File source) {
		file = cacheLookup(dirList, file); 
		if (file == null) return false;
		if (source != null) {
			source = cacheLookup(dirList, source);
			if (!gridss.Defaults.IGNORE_FILE_TIMESTAMPS && source != null && source.exists() && file.lastModified() < source.lastModified()) {
				log.info(source, " has a more recent timestamp than ", file, ". Considering ", file, " out of date.");
				return false;
			}
		}
		return true;
	}
	private static File cacheLookup(Map<String, File[]> dirList, File file) {
		file = file.getAbsoluteFile();
		File dir = file.getParentFile();
		if (dir == null) {
			// bypass cache
			log.debug("Not caching " + file.toString());
			return file;
		}
		String dirKey = dir.getPath();
		File[] dirContents = dirList.get(dirKey);
		if (dirContents == null) {
			dirContents = dir.listFiles();
			dirList.put(dirKey, dirContents);
		}
		if (dirContents != null) { 
			for (File f : dirContents) {
				if (f.getName().equals(file.getName())) {
					return f;
				}
			}
		}
		return null;
	}
}
