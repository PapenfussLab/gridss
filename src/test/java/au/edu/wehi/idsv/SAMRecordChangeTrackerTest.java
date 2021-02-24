package au.edu.wehi.idsv;

import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMRecord;
import org.junit.Assert;
import org.junit.Test;

public class SAMRecordChangeTrackerTest extends TestHelper {
    @Test
    public void should_track_changes() {
        SAMRecord r1 = Read(0, 1, "10M");
        SAMRecord r2 = Read(0, 2, "10M");
        SAMRecordChangeTracker ct = new SAMRecordChangeTracker();
        SAMRecordChangeTracker.TrackedFragment tf = ct.startTrackedChanges(ImmutableList.of(r1));
        r1.setAttribute("zz", 1);
        r1.setCigarString("5M5S");
        ct.processTrackedChanges(tf, ImmutableList.of(r1, r2));
        SAMRecordChangeTracker.Changes c = ct.getChanges();
        Assert.assertEquals(2, c.totalAlignments);
        Assert.assertEquals(1, c.addedAlignments);
        Assert.assertEquals(1, c.updatedCigar);
        Assert.assertEquals(1, c.updatedTag.getLong("zz"));
    }
}