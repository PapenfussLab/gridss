package au.edu.wehi.idsv.sam;

import au.edu.wehi.idsv.util.GroupingIterator;
import au.edu.wehi.idsv.util.MessageThrottler;
import com.google.common.collect.Iterators;
import com.google.common.collect.Ordering;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;

import java.util.Iterator;
import java.util.List;
import java.util.Set;

/**
 * Ensures that every aligned record has an associated NM tag 
 * 
 * @author Daniel Cameron
 *
 */
public class TemplateTagsIterator implements Iterator<List<SAMRecord>> {
	private static final Log log = Log.getInstance(TemplateTagsIterator.class);
	private final Set<String> tags;
	private final Iterator<List<SAMRecord>> it;
	private final boolean softenHardClips;
	private final boolean fixMates;
	private final boolean fixDuplicates;
	private final boolean fixSA;
	private final boolean fixTruncated;
	private final boolean recalculateSupplementary;

	public TemplateTagsIterator(Iterator<List<SAMRecord>> it, boolean softenHardClips, boolean fixMates, boolean fixDuplicates, boolean fixSA, boolean fixTruncated, boolean recalculateSupplementary, Set<String> tags) {
		this.it = Iterators.peekingIterator(it);
		this.softenHardClips = softenHardClips;
		this.fixMates = fixMates;
		this.fixDuplicates = fixDuplicates;
		this.fixSA = fixSA;
		this.fixTruncated = fixTruncated;
		this.recalculateSupplementary = recalculateSupplementary;
		this.tags = tags;
	}

	public static Iterator<List<SAMRecord>> withGrouping(Iterator<SAMRecord> it) {
		return new GroupingIterator(it, Ordering.natural().onResultOf((SAMRecord r) -> r.getReadName()));
	}

	@Override
	public boolean hasNext() {
		return it.hasNext();
	}

	@Override
	public List<SAMRecord> next() {
		List<SAMRecord> records = it.next();
		if (records.size() > 16) {
			if (!MessageThrottler.Current.shouldSupress(log, "many read alignments")) {
				log.warn(String.format("Found %d records with read name \"%s\". GRIDSS requires read names be unique", records.size(), records.get(0).getReadName()));
			}
		}
		SAMRecordUtil.calculateTemplateTags(records, tags, softenHardClips, fixMates, fixDuplicates, fixSA, fixTruncated, recalculateSupplementary);
		return records;
	}
}
