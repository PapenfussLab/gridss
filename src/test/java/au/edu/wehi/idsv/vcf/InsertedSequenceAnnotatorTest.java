package au.edu.wehi.idsv.vcf;

import au.edu.wehi.idsv.TestHelper;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMRecord;
import org.junit.Assert;
import org.junit.Test;

import java.util.List;

public class InsertedSequenceAnnotatorTest extends TestHelper {
    @Test
    public void should_replace_pipes_with_underscores() {
        SAMRecord r = new SAMRecord(null);
        r.setReferenceName("with|pipe");
        r.setAlignmentStart(1);
        r.setCigarString("50M");
        r.setMappingQuality(5);
        r.setAttribute("XA", "chr|piped,-5,50M,6;");
        List<String> result = InsertedSequenceAnnotator.writeAlignmentAnnotation(ImmutableList.of(r));
        Assert.assertEquals(2, result.size());
        Assert.assertEquals("with_pipe:1|+|50M|5", result.get(0));
        Assert.assertEquals("chr_piped:5|-|50M|", result.get(1));
    }
}