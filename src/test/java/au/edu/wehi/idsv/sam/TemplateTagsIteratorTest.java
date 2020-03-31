package au.edu.wehi.idsv.sam;

import au.edu.wehi.idsv.TestHelper;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import org.junit.Test;

import java.io.File;
import java.util.List;

import static org.junit.Assert.assertEquals;

public class TemplateTagsIteratorTest extends TestHelper {
    @Test
    public void should_jointly_process_tags() {
        List<SAMRecord> rec = getRecords(new File("src/test/resources/flag.sam"));
        TemplateTagsIterator tti = new TemplateTagsIterator(rec.iterator(), true, true, true, true, true, true, Sets.newHashSet(
                SAMTag.NM.name(),
                SAMTag.SA.name(),
                SAMTag.R2.name(),
                SAMTag.MC.name(),
                SAMTag.MQ.name()));
        List<SAMRecord> out = Lists.newArrayList(tti);
        assertEquals(81, out.get(0).getFlags());
        assertEquals(161, out.get(1).getFlags());
        assertEquals(2209, out.get(2).getFlags());
        assertEquals("110S40M", out.get(2).getCigarString());
    }

}