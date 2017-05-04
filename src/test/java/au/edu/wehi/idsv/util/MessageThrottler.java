package au.edu.wehi.idsv.util;

import java.util.HashMap;

public class MessageThrottler {
	public static final MessageThrottler Current = new MessageThrottler(gridss.Defaults.SUPPRESS_DATA_ERROR_MESSAGES_AFTER);
	private HashMap<String, Integer> counts = new HashMap<>();
	private final int threshold;
	public MessageThrottler(int suppressDataErrorMessagesAfter) {
		this.threshold = suppressDataErrorMessagesAfter;
	}
	public boolean shouldSupress(htsjdk.samtools.util.Log log, String messageName) {
		Integer countObj = counts.get(messageName);
		int count = countObj == null ? 0 : countObj.intValue();
		count++;
		counts.put(messageName, count);
		if (count == threshold) {
			log.info(String.format("Supressing messages regarding %s", messageName));
		}
		return count < threshold;
	}
}
