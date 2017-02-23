package gridss;

public class Defaults {
	public static final boolean IGNORE_FILE_TIMESTAMPS;
	public static final int ASYNC_BUFFERS;
	public static final int ASYNC_BUFFER_SIZE;
	public static final boolean DELETE_TEMPORARY_FILES;
	public static final int SUPPRESS_DATA_ERROR_MESSAGES_AFTER;
	public static final boolean WRITE_ZERO_OR_EMTPY_VCF_FIELDS;
	static {
		IGNORE_FILE_TIMESTAMPS = Boolean.valueOf(System.getProperty("gridss.ignoreTimestamps", "true"));
		ASYNC_BUFFERS = Integer.parseInt(System.getProperty("gridss.async.buffers", "2"));
		ASYNC_BUFFER_SIZE = Integer.parseInt(System.getProperty("gridss.async.buffers", "300"));
		DELETE_TEMPORARY_FILES = !Boolean.valueOf(System.getProperty("gridss.keepTempFiles", "false"));
		SUPPRESS_DATA_ERROR_MESSAGES_AFTER = Integer.parseInt(System.getProperty("gridss.logSpamLimit", "100"));
		WRITE_ZERO_OR_EMTPY_VCF_FIELDS = Boolean.valueOf(System.getProperty("gridss.writeZeroOrEmptyVcfFields", "true"));
	}
}
