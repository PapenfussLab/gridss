package au.edu.wehi.idsv.repeatmasker;

import au.edu.wehi.idsv.TestHelper;
import com.google.common.collect.Lists;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.readers.LineIterator;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;

public class RepeatMaskerCodecTest extends TestHelper {
    @Test
    public void should_parse_summary_long_query_sequence() throws IOException {
        File file = new File("src/test/resources/repeatmasker/rmtest.fa.out");
        try (AbstractFeatureReader<RepeatMaskerFeature, LineIterator> reader = AbstractFeatureReader.getFeatureReader(file.getPath(), new RepeatMaskerCodec(), false)) {
            Assert.assertNull(reader.getHeader());
            List<RepeatMaskerFeature> features = Lists.newArrayList((Iterable<? extends RepeatMaskerFeature>) reader.iterator());
            Assert.assertEquals(2, features.size());
            RepeatMaskerFeature f = features.get(0);
            Assert.assertEquals(327, f.getSmithWatermanScore());
            Assert.assertEquals(0.25, f.getRepeatAlignmentSummaryInformation().getPercentageSubstituted(), 0);
            Assert.assertEquals(0.5, f.getRepeatAlignmentSummaryInformation().getPercentageDeleted(), 0);
            Assert.assertEquals(0.75, f.getRepeatAlignmentSummaryInformation().getPercentageInserted(), 0);
            Assert.assertEquals("AB007990REALLYLONGCONTIGNAME_ASDF_ASDFSDF_ASD_FA_S", f.getContig());
            Assert.assertEquals(4, f.getStart());
            Assert.assertEquals(357, f.getEnd());
            Assert.assertEquals(3, f.getRepeatAlignmentSummaryInformation().getBasesInQueryPastMatch());
            Assert.assertEquals(Strand.POSITIVE, f.getStrand());
            Assert.assertEquals("(A)n", f.getRepeatType());
            Assert.assertEquals("Simple_repeat", f.getRepeatClass());
            Assert.assertEquals(2, f.getRepeatAlignmentSummaryInformation().getMatchStart());
            Assert.assertEquals(354, f.getRepeatAlignmentSummaryInformation().getMatchEnd());
            Assert.assertEquals(0, f.getRepeatAlignmentSummaryInformation().getBasesInRepeatPastMatch());
            Assert.assertEquals("1", f.getUniqueID());
            Assert.assertEquals("3S354M", f.getRepeatAlignmentInformation(true).getCigar().toString());
            Assert.assertEquals(2, f.getRepeatAlignmentInformation(true).getRepeatStart());
            Assert.assertEquals(0, f.getRepeatAlignmentInformation(true).getNestedBases());
            Assert.assertNull(f.getRepeatAlignmentInformation(false));

            // complement sequence reverses order (start,end,offset) tuple
            f = features.get(1);
            Assert.assertEquals(208, f.getRepeatAlignmentSummaryInformation().getMatchStart());
            Assert.assertEquals(537, f.getRepeatAlignmentSummaryInformation().getMatchEnd());
            Assert.assertEquals(12, f.getRepeatAlignmentSummaryInformation().getBasesInRepeatPastMatch());
        }
    }
    @Test
    public void should_parse_alignment_file() throws IOException {
        File file = new File("src/test/resources/repeatmasker/rmtest.fa.cat");
        try (AbstractFeatureReader<RepeatMaskerFeature, LineIterator> reader = AbstractFeatureReader.getFeatureReader(file.getPath(), new RepeatMaskerCodec(), false)) {
            Assert.assertNull(reader.getHeader());
            List<RepeatMaskerFeature> features = Lists.newArrayList((Iterable<? extends RepeatMaskerFeature>) reader.iterator());

            RepeatMaskerFeature f = features.get(0);
            Assert.assertEquals("(TCTAA)n", f.getRepeatType());
            Assert.assertEquals("Simple_repeat", f.getRepeatClass());
            Assert.assertEquals("m_b1s252i11", f.getUniqueID());
            Assert.assertEquals(1, f.getRepeatAlignmentInformation(true).getRepeatStart());
            //  AC_000192             33 TCTACTCTAAAACTCTTGTAGTTTAAATCTAATCTAATCTAAACT 77
            //                               v   -  v  vv v  i i -                 v
            //  (TCTAA)n#Simp          1 TCTAATCT-AATCTAATCTAATCT-AATCTAATCTAATCTAATCT 43
            //                           ====x===i==x==xx=x==x=x=i=================x==
            //                           4   13  12 12 2 112 1111117               12
            //
            Assert.assertEquals("32S4=1X3=1I2=1X2=2X1=1X2=1X1=1X1=1I17=1X2=", f.getRepeatAlignmentInformation(true).getCigar().toString());
            Assert.assertEquals(0, f.getRepeatAlignmentInformation(true).getNestedBases());

            // multi-line alignment
            f = features.get(12);
            Assert.assertEquals("AY849321", f.getContig());
            Assert.assertEquals("53S5=1I4=1I2=1X1=1X6=1I4=1X6=5X4=1X3=1D4=1X6=1I1X2=1X4=", f.getRepeatAlignmentInformation(true).getCigar().toString());
            Assert.assertEquals(1, f.getRepeatAlignmentInformation(true).getRepeatStart());
            //  AY849321              54 TTTTCTTTTTACTGTACTTTTCTTTTTGTTTTCTACATGTTTCATTT-TT 102
            //                                -    -  v v      -    v      viviv    v   -
            //  (TTTTC)n#Simp          1 TTTTC-TTTT-CTTTTCTTTTC-TTTTCTTTTCTTTTCTTTTCTTTTCTT 47
            //
            //  AY849321             103 TTGTTTTCTGGTTATTTT 120
            //                             v      -v  v
            //  (TTTTC)n#Simp         48 TTCTTTTCT-TTTCTTTT 64

            // no starting soft clip
            f = features.stream().filter(x -> x.getUniqueID().equals("m_b4s252i6")).findFirst().get();
            Assert.assertEquals("27=1X3=1X7=", f.getRepeatAlignmentInformation(true).getCigar().toString());
        }
    }
}