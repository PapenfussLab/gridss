/**
 * 
 */
package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import it.unimi.dsi.fastutil.objects.Object2LongMap;
import it.unimi.dsi.fastutil.objects.Object2LongOpenHashMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;

/**
 * @author Daniel Cameron
 *
 */
public class SAMRecordChangeTracker {
	private static final String TRANSIENT_TRACKING_ATTRIBUTE_NAME = "scrt";
	private static final int N_FLAGS = 16;
	private final Changes changes = new Changes();

	public Changes getChanges() { return changes; }

	public static class Changes {
		public long totalFragments = 0;
		public long totalAlignments = 0;
		public long removedAlignments = 0;
		public long addedAlignments = 0;
		public long updatedQName = 0;
		public final long[] setFlag = new long[N_FLAGS];
		public final long[] clearedFlag = new long[N_FLAGS];
		public final long[] updatedFlag = new long[N_FLAGS];
		public long updatedRName = 0;
		public long updatedPos = 0;
		public long updatedMapq = 0;
		public long updatedCigar = 0;
		public long updatedRNext = 0;
		public long updatedPNext = 0;
		public long updatedTLen = 0;
		public long updatedSeq = 0;
		public long updatedQual = 0;
		public final Object2LongOpenHashMap<String> updatedTag = new Object2LongOpenHashMap<>();
		public void accumulate(Changes c) {
			totalFragments += c.totalFragments;
			totalAlignments += c.totalAlignments;
			removedAlignments += c.removedAlignments;
			addedAlignments += c.addedAlignments;
			updatedQName += c.updatedQName;
			updatedRName += c.updatedRName;
			updatedPos += c.updatedPos;
			updatedMapq += c.updatedMapq;
			updatedCigar += c.updatedCigar;
			updatedRNext += c.updatedRNext;
			updatedPNext += c.updatedPNext;
			updatedTLen += c.updatedTLen;
			updatedSeq += c.updatedSeq;
			updatedQual += c.updatedQual;
			for (int i = 0; i < N_FLAGS; i++) {
				setFlag[i] += c.setFlag[i];
				clearedFlag[i] += c.clearedFlag[i];
				updatedFlag[i] += c.updatedFlag[i];
			}
			for (Object2LongMap.Entry<String> tag : c.updatedTag.object2LongEntrySet()) {
				updatedTag.addTo(tag.getKey(), tag.getLongValue());
			}
		}
		public boolean hasChanges() {
			if (removedAlignments > 0) return true;
			if (addedAlignments > 0) return true;
			if (updatedQName > 0) return true;
			for (int i = 0; i < N_FLAGS; i++) {
				if (setFlag[i] > 0) return true;
				if (clearedFlag[i] > 0) return true;
				if (updatedFlag[i] > 0) return true;
			}
			if (updatedRName > 0) return true;
			if (updatedPos > 0) return true;
			if (updatedMapq > 0) return true;
			if (updatedCigar > 0) return true;
			if (updatedRNext > 0) return true;
			if (updatedPNext > 0) return true;
			if (updatedTLen > 0) return true;
			if (updatedSeq > 0) return true;
			if (updatedQual > 0) return true;
			if (updatedTag.size() > 0) return true;
			return false;
		}
	}
	public static class TrackedFragment {
		private List<SAMRecord> before;
		public TrackedFragment(List<SAMRecord> records) {
			try {
				this.before = new ArrayList<>(records.size());
				int i = 0;
				for (SAMRecord r : records) {
					r.setTransientAttribute(TRANSIENT_TRACKING_ATTRIBUTE_NAME, i);
					SAMRecord mirror = (SAMRecord)r.clone();
					before.add(mirror);
				}
			} catch (CloneNotSupportedException e) {
			}
		}
	}
	public TrackedFragment startTrackedChanges(List<SAMRecord> before) {
		return new TrackedFragment(before);
	}
	public void processTrackedChanges(TrackedFragment before, List<SAMRecord> after) {
		Changes c = new Changes();
		for (int i = 0; i < before.before.size(); i++) {
			SAMRecord b = before.before.get(i);
			boolean matchFound = false;
			for (SAMRecord a : after) {
				Object afterKey = a.getTransientAttribute(TRANSIENT_TRACKING_ATTRIBUTE_NAME);
				if (afterKey instanceof Integer && (int)afterKey == i) {
					matchFound= true;
					processAlignment(c, b, a);
				}
			}
			if (!matchFound) {
				processAlignment(c, b, null);
			}
		}
		for (SAMRecord r : after) {
			if (r.getTransientAttribute(TRANSIENT_TRACKING_ATTRIBUTE_NAME) == null) {
				processAlignment(c, null, r);
			}
		}
		c.totalFragments += 1;
		// Minimise the time spent in the lock by only doing a full update when there are actually changed records
		if (c.hasChanges()) {
			synchronized (this.changes) {
				this.changes.accumulate(c);
			}
		} else {
			synchronized (this.changes) {
				this.changes.totalFragments += c.totalFragments;
				this.changes.totalAlignments += c.totalAlignments;
			}
		}
    }

	private static void processAlignment(Changes c, SAMRecord b, SAMRecord a) {
		if (b == null && a == null) return;
		c.totalAlignments++;
		if (b == null) {
			c.addedAlignments++;
		} else if (a == null) {
			c.removedAlignments++;
		} else {
			if (!Objects.equals(a.getReadName(), b.getReadName())) {
				c.updatedQName++;
			}
			int aFlag = a.getFlags();
			int bFlag = b.getFlags();
			for (int i = 0; i < N_FLAGS; i++) {
				if ((aFlag & (1 << i)) != (bFlag & (1 << i))) {
					c.updatedFlag[i]++;
					if ((aFlag & (1 << i)) == 0) {
						c.clearedFlag[i]++;
					} else {
						c.setFlag[i]++;
					}
				}
			}
			if (a.getReferenceIndex() != b.getReferenceIndex()) {
				c.updatedRName++;
			}
			if (a.getAlignmentStart() != b.getAlignmentStart()) {
				c.updatedPos++;
			}
			if (a.getMappingQuality() != b.getMappingQuality()) {
				c.updatedPos++;
			}
			if (!Objects.equals(a.getMappingQuality(), b.getMappingQuality())) {
				c.updatedMapq++;
			}
			if (!Objects.equals(a.getCigar(), b.getCigar())) {
				c.updatedCigar++;
			}
			if (a.getMateReferenceIndex() != b.getMateReferenceIndex()) {
				c.updatedRNext++;
			}
			if (a.getMateAlignmentStart() != b.getMateAlignmentStart()) {
				c.updatedPNext++;
			}
			if (a.getInferredInsertSize() != b.getInferredInsertSize()) {
				c.updatedTLen++;
			}
			if (!Arrays.equals(a.getReadBases(), b.getReadBases())) {
				c.updatedSeq++;
			}
			if (!Arrays.equals(a.getBaseQualities(), b.getBaseQualities())) {
				c.updatedQual++;
			}
			for (SAMRecord.SAMTagAndValue attr : a.getAttributes()) {
				Object bValue = b.getAttribute(attr.tag);
				if (bValue == null) {
					// tag added
					c.updatedTag.addTo(attr.tag, 1);
				} else if (!areAttributesEqual(attr.value, b.getAttribute(attr.tag))) {
					// tag updated
					c.updatedTag.addTo(attr.tag, 1);
				}
			}
			for (SAMRecord.SAMTagAndValue attr : b.getAttributes()) {
				if (!a.hasAttribute(attr.tag)) {
					// tag removed
					c.updatedTag.addTo(attr.tag, 1);
				}
			}
		}
	}
	private static boolean areAttributesEqual(Object a, Object b) {
		if (a == null && b == null) return true;
		if (a == null || b == null) return false;
		// The HTSJDK API isn't good with consistent type persistence
		// (e.g. ints being returned as strings)
		// so we just hackily convert everything to string instead
		// of actually dealing with all the underlying types
		return a.toString().equals(b.toString());
	}
	public void writeSummary(File output) throws IOException {
		List<String> lines = new ArrayList<>();
		lines.add("totalFragments,"+ changes.totalFragments);
		lines.add("totalAlignments,"+ changes.totalAlignments);
		lines.add("removedAlignments,"+ changes.removedAlignments);
		lines.add("addedAlignments,"+ changes.addedAlignments);
		lines.add("updatedQName,"+ changes.updatedQName);
		for (int i = 0; i < N_FLAGS; i++) {
			lines.add("setFlag_"+ (1 << i) +","+ changes.setFlag[i]);
			lines.add("clearedFlag_"+ (1 << i) +","+ changes.clearedFlag[i]);
			lines.add("updatedFlag_"+ (1 << i) +","+ changes.updatedFlag[i]);
		}
		lines.add("updatedRName,"+ changes.updatedRName);
		lines.add("updatedPos,"+ changes.updatedPos);
		lines.add("updatedMapq,"+ changes.updatedMapq);
		lines.add("updatedCigar,"+ changes.updatedCigar);
		lines.add("updatedRNext,"+ changes.updatedRNext);
		lines.add("updatedPNext,"+ changes.updatedPNext);
		lines.add("updatedTLen,"+ changes.updatedTLen);
		lines.add("updatedSeq,"+ changes.updatedSeq);
		lines.add("updatedQual,"+ changes.updatedQual);
		for (Object2LongMap.Entry<String> entry : changes.updatedTag.object2LongEntrySet()) {
			lines.add("updatedTag_"+entry.getKey()+","+entry.getLongValue());
		}
		Files.write(output.toPath(), lines);
	}
}

