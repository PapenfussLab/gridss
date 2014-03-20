package au.edu.wehi.socrates;

/**
 * Set of reads providing support for a given position
 * @author Daniel Cameron
 *
 */
public class SAMRecordWorkingSet {
	private int scIntervalOffsetStart = 1;
	private int scIntervalOffsetEnd = -1;
	/**
	 * Inver
	 * @author Daniel Cameron
	 *
	 */
	public class EvidenceSupportInterval {}
	/**
	 * Adds the given record to the working set
	 * @param record
	 * @return
	 */
	public Integer addToWorkingSet(SAMRecord record);
	public int getWorkingSetStartPosition(SAMRecord record);
	public int getWorkingSetEndPosition(SAMRecord record);
	public void setPosition(int position);
	public int getPosition();
}
