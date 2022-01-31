package gridss;

import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.IntermediateFilesTest;
import au.edu.wehi.idsv.NonReferenceReadPair;
import au.edu.wehi.idsv.SoftClipEvidence;
import au.edu.wehi.idsv.sam.SAMRecordUtil;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import gridss.cmdline.CommandLineProgramHelper;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.ValidationStringency;
import org.junit.Assert;
import org.junit.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import static org.junit.Assert.assertEquals;

public class ComputeSamTagsTest extends IntermediateFilesTest {
    @Test
    public void should_update_NM_R2() throws IOException {
        List<DirectedEvidence> in = new ArrayList<DirectedEvidence>();
        in.add(SCE(FWD, withSequence("AACCGGTTCTA", Read(0, 15, "5M6S"))));
        in.add(SCE(FWD, withSequence("AACCGGTTCTA", Read(0, 15, "6M5S"))));
        in.add(NRRP(withSequence("AACCGGTTCTA", DP(0, 1, "11M", true, 1, 100, "11M", true))));
        in.add(NRRP(withSequence("AACCGGTTCTA", DP(0, 2, "11M", true, 1, 100, "11M", false))));

        List<SAMRecord> insam = in.stream().flatMap(de -> {
            if (de instanceof SoftClipEvidence) return ImmutableList.<SAMRecord>of(((SoftClipEvidence)de).getSAMRecord()).stream();
            if (de instanceof NonReferenceReadPair) return ImmutableList.<SAMRecord>of(
                    ((NonReferenceReadPair)de).getNonReferenceRead(),
                    ((NonReferenceReadPair)de).getLocalledMappedRead()).stream();
            throw new RuntimeException("NYI");
        }).collect(Collectors.toList());
        insam.addAll(Lists.newArrayList(RP(0, 1, 100, 10)));
        insam.forEach(r -> r.setAttribute("NM", null));
        Collections.sort(insam, new SAMRecordQueryNameComparator());
        createBAM(input, SAMFileHeader.SortOrder.queryname, insam);
        CommandLineProgramHelper cmd = new CommandLineProgramHelper(new ComputeSamTags());
        cmd.addArg("INPUT", input);
        cmd.addArg("OUTPUT" , output);
        cmd.addArg("REFERENCE_SEQUENCE" , reference);
        cmd.addArg("TAGS" , "null");
        cmd.addArg("TAGS" , "R2");
        cmd.addArg("TAGS" , "NM");
        cmd.addArg("ASSUME_SORTED" , true);
        Assert.assertEquals(0, cmd.run());
        List<SAMRecord> out = getRecords(output);
        Assert.assertEquals(8, (int)out.get(0).getIntegerAttribute("NM"));
        Assert.assertEquals("AACCGGTTCTA", out.get(0).getStringAttribute("R2"));

    }
    @Test
    public void should_fix_unmapped_primary_alignment() throws IOException {
        List<DirectedEvidence> in = new ArrayList<DirectedEvidence>();
        for (boolean withSA : new boolean[] { true, false}) {
            SAMRecord[] reads = withName("r", withSequence("AACCGGTTCTA",
                    Read(0, 15, "5M6S"),
                    Read(1, 15, "5S6M")));
            if (withSA) {
                reads = withAttr("SA", "polyA,100,+,2S8M,0,0", reads);
            }
            reads[1].setSupplementaryAlignmentFlag(true);
            reads[0].setReadUnmappedFlag(true);

            List<SAMRecord> insam = Lists.newArrayList(reads);
            Collections.sort(insam, new SAMRecordQueryNameComparator());
            createBAM(input, SAMFileHeader.SortOrder.queryname, insam);
            CommandLineProgramHelper cmd = new CommandLineProgramHelper(new ComputeSamTags());
            cmd.addArg("INPUT", input);
            cmd.addArg("OUTPUT", output);
            cmd.addArg("REFERENCE_SEQUENCE", reference);
            cmd.addArg("ASSUME_SORTED", true);
            cmd.getProgram().VALIDATION_STRINGENCY = ValidationStringency.LENIENT;
            Assert.assertEquals(0, cmd.run());
            List<SAMRecord> out = getRecords(output);
            Assert.assertEquals(1, out.size());
            Assert.assertEquals((Integer)1, out.get(0).getReferenceIndex());
            Assert.assertFalse(out.get(0).getSupplementaryAlignmentFlag());
        }
    }
    @Test
    public void should_remove_records_when_necessary() {
        List<SAMRecord> list = new ArrayList<>();
        list.add(Read(0, 1, "10S90M"));
        list.add(Read(0, 100, "90M10S"));
        list.forEach(r -> r.setReadName("r"));
        for (int i = 1; i < list.size(); i++) {
            list.get(i).setSupplementaryAlignmentFlag(true);
        }
        Collections.sort(list, new SAMRecordQueryNameComparator());
        createBAM(input, SAMFileHeader.SortOrder.queryname, list);
        CommandLineProgramHelper cmd = new CommandLineProgramHelper(new ComputeSamTags());
        cmd.addArg("INPUT", input);
        cmd.addArg("OUTPUT", output);
        cmd.addArg("REFERENCE_SEQUENCE", reference);
        cmd.addArg("ASSUME_SORTED", true);
        cmd.getProgram().VALIDATION_STRINGENCY = ValidationStringency.LENIENT;
        Assert.assertEquals(0, cmd.run());
        List<SAMRecord> out = getRecords(output);
        Assert.assertEquals(1, out.size());
    }
}