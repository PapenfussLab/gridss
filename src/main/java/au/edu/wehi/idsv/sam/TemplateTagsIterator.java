package au.edu.wehi.idsv.sam;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;

import java.util.ArrayDeque;
import java.util.Iterator;
import java.util.Queue;

import com.google.common.collect.Iterators;
import com.google.common.collect.PeekingIterator;

/**
 * Ensures that every aligned record has an associated NM tag 
 * 
 * @author Daniel Cameron
 *
 */
public class TemplateTagsIterator implements Iterator<SAMRecord> {
	private final PeekingIterator<SAMRecord> it;
	private Queue<SAMRecord> queue = new ArrayDeque<>();
	public TemplateTagsIterator(Iterator<SAMRecord> it) {
		this.it = Iterators.peekingIterator(it);
	}
	private void ensureQueue() {
		List<SAMRecord>
		if (it.hasNext()) {
			String readname = it.peek().getReadName();
			while (readname != null && readname.equals(it.peek().getReadName())) {
				
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
