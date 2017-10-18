package gridss;

public class Defaults {
	public static final boolean IGNORE_FILE_TIMESTAMPS;
	public static final int ASYNC_BUFFERS;
	public static final int ASYNC_BUFFER_SIZE;
	public static final boolean DELETE_TEMPORARY_FILES;
	public static final int SUPPRESS_DATA_ERROR_MESSAGES_AFTER;
	public static final boolean WRITE_ZERO_OR_EMTPY_VCF_FIELDS;
	/**
	 * Output to temporary file. This allows restarting when a process is killed
	 * without having to worry about partially written output files being used
	 * when the process is restarted. 
	 */
	public static final boolean OUTPUT_TO_TEMP_FILE;
	/**
	 * Perform defensive garbage collection to reduce the chance of running out of file handles.
	 * See http://stackoverflow.com/questions/2972986/how-to-unmap-a-file-from-memory-mapped-using-filechannel-in-java
	 */
	public static final boolean DEFENSIVE_GC;
	static {
		IGNORE_FILE_TIMESTAMPS = Boolean.valueOf(System.getProperty("gridss.ignoreTimestamps", "true"));
		ASYNC_BUFFERS = Integer.parseInt(System.getProperty("gridss.async.buffers", "2"));
		ASYNC_BUFFER_SIZE = Integer.parseInt(System.getProperty("gridss.async.buffers", "300"));
		DELETE_TEMPORARY_FILES = !Boolean.valueOf(System.getProperty("gridss.keepTempFiles", "false"));
		SUPPRESS_DATA_ERROR_MESSAGES_AFTER = Integer.parseInt(System.getProperty("gridss.logSpamLimit", "100"));
		WRITE_ZERO_OR_EMTPY_VCF_FIELDS = Boolean.valueOf(System.getProperty("gridss.writeZeroOrEmptyVcfFields", "true"));
		DEFENSIVE_GC = Boolean.valueOf(System.getProperty("gridss.defensiveGC", "false"));
		OUTPUT_TO_TEMP_FILE = Boolean.valueOf(System.getProperty("gridss.output_to_temp_file", "false"));
	}
}
