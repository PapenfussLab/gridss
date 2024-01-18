package gridss.cmdline;

import au.edu.wehi.idsv.IntermediateFilesTest;
import au.edu.wehi.idsv.ProcessingContext;
import com.google.common.collect.ImmutableList;
import gridss.SanityCheckEvidence;
import htsjdk.samtools.SAMFileHeader;
import org.junit.Test;

import java.io.File;

public class FullEvidenceCommandLineProgramTest extends IntermediateFilesTest {
    private FullEvidenceCommandLineProgram create(ProcessingContext pc) {
        SanityCheckEvidence cmd = new SanityCheckEvidence();
        cmd.WORKING_DIR = testFolder.getRoot();
        cmd.setReference(pc.getReferenceFile());
        cmd.setReference(pc.getReference());
        cmd.setContext(pc);
        return cmd;
    }
    @Test
    public void should_be_happy_when_all_labels_assembled(){
        ProcessingContext pc = getCommandlineContext();
        File n = new File(testFolder.getRoot(), "n.bam");
        File t = new File(testFolder.getRoot(), "t.bam");
        File a1 = new File(testFolder.getRoot(), "a1.bam");
        File a2 = new File(testFolder.getRoot(), "a2.bam");
        createBAM(n, SAMFileHeader.SortOrder.coordinate);
        createBAM(t, SAMFileHeader.SortOrder.coordinate);
        SAMFileHeader h1 = AES().getHeader();
        SAMFileHeader h2 = AES().getHeader();
        h1.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        h1.setComments(ImmutableList.of("gridss_input_category=Normal"));
        h2.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        h2.setComments(ImmutableList.of("gridss_input_category=Tumour"));
        createBAM(a1, h1);
        createBAM(a2, h2);
        FullEvidenceCommandLineProgram cmd = create(pc);
        cmd.INPUT = ImmutableList.of(n, t);
        cmd.INPUT_LABEL = ImmutableList.of("Normal", "Tumour");
        cmd.ASSEMBLY = ImmutableList.of(a1, a2);
        cmd.SAMPLE_NAMES = ImmutableList.of();
        cmd.doWork();
    }
    @Test(expected=IllegalArgumentException.class)
    public void should_complain_about_missing_assembly_label(){
        ProcessingContext pc = getCommandlineContext();
        File n = new File(testFolder.getRoot(), "n.bam");
        File t = new File(testFolder.getRoot(), "t.bam");
        File a1 = new File(testFolder.getRoot(), "a1.bam");
        File a2 = new File(testFolder.getRoot(), "a2.bam");
        createBAM(n, SAMFileHeader.SortOrder.coordinate);
        createBAM(t, SAMFileHeader.SortOrder.coordinate);
        SAMFileHeader h1 = AES().getHeader();
        SAMFileHeader h2 = AES().getHeader();
        h1.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        h1.setComments(ImmutableList.of("gridss_input_category=Normal"));
        h2.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        h2.setComments(ImmutableList.of("gridss_input_category=Normal"));
        createBAM(a1, h1);
        createBAM(a2, h2);
        FullEvidenceCommandLineProgram cmd = create(pc);
        cmd.INPUT = ImmutableList.of(n, t);
        cmd.INPUT_LABEL = ImmutableList.of("Normal", "Tumour");
        cmd.ASSEMBLY = ImmutableList.of(a1, a2);
        cmd.SAMPLE_NAMES = ImmutableList.of();
        cmd.doWork();
    }
    @Test(expected=IllegalArgumentException.class)
    public void should_complain_about_duplicate_assembly_label() {
        ProcessingContext pc = getCommandlineContext();
        File n = new File(testFolder.getRoot(), "n.bam");
        File t = new File(testFolder.getRoot(), "t.bam");
        File a1 = new File(testFolder.getRoot(), "a1.bam");
        File a2 = new File(testFolder.getRoot(), "a2.bam");
        createBAM(n, SAMFileHeader.SortOrder.coordinate);
        createBAM(t, SAMFileHeader.SortOrder.coordinate);
        SAMFileHeader h1 = AES().getHeader();
        SAMFileHeader h2 = AES().getHeader();
        h1.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        h1.setComments(ImmutableList.of("gridss_input_category=Normal"));
        h2.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        h2.setComments(ImmutableList.of("gridss_input_category=Normal"));
        createBAM(a1, h1);
        createBAM(a2, h2);
        FullEvidenceCommandLineProgram cmd = create(pc);
        cmd.INPUT = ImmutableList.of(n, t);
        cmd.INPUT_LABEL = ImmutableList.of("Normal", "Tumour");
        cmd.ASSEMBLY = ImmutableList.of(a1, a2);
        cmd.SAMPLE_NAMES = ImmutableList.of();
        cmd.doWork();
    }
    @Test(expected=IllegalArgumentException.class)
    public void should_complain_about_duplicate_assembly_name() {
        ProcessingContext pc = getCommandlineContext();
        File n = new File(testFolder.getRoot(), "n.bam");
        File t = new File(testFolder.getRoot(), "t.bam");
        File a1 = new File(testFolder.getRoot(), "a1.bam");
        File a2 = new File(testFolder.getRoot(), "a1.bam");
        createBAM(n, SAMFileHeader.SortOrder.coordinate);
        createBAM(t, SAMFileHeader.SortOrder.coordinate);
        SAMFileHeader h1 = AES().getHeader();
        SAMFileHeader h2 = AES().getHeader();
        h1.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        h1.setComments(ImmutableList.of("gridss_input_category=Normal"));
        h2.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        h2.setComments(ImmutableList.of("gridss_input_category=Tumour"));
        createBAM(a1, h1);
        createBAM(a2, h2);
        FullEvidenceCommandLineProgram cmd = create(pc);
        cmd.INPUT = ImmutableList.of(n, t);
        cmd.INPUT_LABEL = ImmutableList.of("Normal", "Tumour");
        cmd.ASSEMBLY = ImmutableList.of(a1, a2);
        cmd.SAMPLE_NAMES = ImmutableList.of();
        cmd.doWork();
    }
}