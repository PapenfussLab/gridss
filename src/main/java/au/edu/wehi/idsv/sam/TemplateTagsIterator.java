package au.edu.wehi.idsv.sam;

import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;
import htsjdk.samtools.SAMRecord;

import java.util.*;

/**
 * Ensures that every aligned record has an associated NM tag 
 * 
 * @author Daniel Cameron
 *
 */
public class TemplateTagsIterator implements Iterator<SAMRecord> {
	private final Set<String> tags;
	private final PeekingIterator<SAMRecord> it;
	private final boolean softenHardClips;
	private final boolean fixMates;
	private final boolean fixDuplicates;
	private final boolean recalculateSupplementary;
	private final Queue<SAMRecord> queue = new ArrayDeque<>();
	public TemplateTagsIterator(Iterator<SAMRecord> it, boolean softenHardClips, boolean fixMates, boolean fixDuplicates, boolean recalculateSupplementary, Set<String> tags) {
		this.it = Iterators.peekingIterator(it);
		this.softenHardClips = softenHardClips;
		this.fixMates = fixMates;
		this.fixDuplicates = fixDuplicates;
		this.recalculateSupplementary = recalculateSupplementary;
		this.tags = tags;
	}
	private void ensureQueue() {
		if (!queue.isEmpty()) return;
		List<SAMRecord> records = new ArrayList<>();
		if (it.hasNext()) {
			String readname = it.peek().getReadName();
			if (readname == null) {
				queue.add(it.next());
			} else {
				while (readname != null && it.hasNext() && readname.equals(it.peek().getReadName())) {
					records.add(it.next());
				}
				SAMRecordUtil.calculateTemplateTags(records, tags, softenHardClips, fixMates, fixDuplicates, recalculateSupplementary);
				queue.addAll(records);
			}
		}
	}
	@Override
	public boolean hasNext() {
		ensureQueue();
		return !queue.isEmpty();
	}
	
	@Override
	public SAMRecord next() {
		ensureQueue();
		return queue.poll();
	}

}
