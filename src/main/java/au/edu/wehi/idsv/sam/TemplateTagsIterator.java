package au.edu.wehi.idsv.sam;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Queue;
import java.util.Set;

import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

/**
 * Ensures that every aligned record has an associated NM tag 
 * 
 * @author Daniel Cameron
 *
 */
public class TemplateTagsIterator implements Iterator<SAMRecord> {
	private final Set<SAMTag> tags;
	private final PeekingIterator<SAMRecord> it;
	private final boolean softenHardClips;
	private final Queue<SAMRecord> queue = new ArrayDeque<>();
	public TemplateTagsIterator(Iterator<SAMRecord> it, boolean softenHardClips, Set<SAMTag> tags) {
		this.it = Iterators.peekingIterator(it);
		this.softenHardClips = softenHardClips;
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
				SAMRecordUtil.calculateTemplateTags(records, tags, softenHardClips);
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
