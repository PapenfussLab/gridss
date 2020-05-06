package gridss;

import au.edu.wehi.idsv.*;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Streams;
import org.apache.commons.math3.util.Pair;
import org.junit.Test;
import org.junit.experimental.categories.Category;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.stream.Collectors;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class GeneratePonBedpeTest extends IntermediateFilesTest {
    private static final File colo829 = new File("src/test/resources/colo829_1.gridss.somatic.vcf");

    private static ProcessingContext getBroadHg19PC() {
        return new ProcessingContext(getFSContext(), Hg19Tests.findBroadHg19Reference(), null, ImmutableList.of(), getConfig());
    }
    @Test
    @Category(Hg19Tests.class)
    public void should_generate_pon() throws IOException {
        File out_bed = new File(testFolder.getRoot(), "out.bed");
        File out_bedpe = new File(testFolder.getRoot(), "out.bedpe");
        new GeneratePonBedpe().instanceMain(new String[] {
                "INPUT=" + colo829.getAbsolutePath(),
                "OUTPUT_BED=" + out_bed.getAbsolutePath(),
                "OUTPUT_BEDPE=" + out_bedpe.getAbsolutePath(),
                "NORMAL_ORDINAL=1",
                "R=" + Hg19Tests.findBroadHg19Reference().getAbsolutePath()
        });
        assertTrue(out_bed.exists());
        assertTrue(out_bedpe.exists());
        List<IdsvVariantContext> vc = getVcf(colo829, getBroadHg19PC(), null);
        List<String> bedpe = Files.readAllLines(out_bedpe.toPath()).stream().filter(s -> !s.startsWith("#")).collect(Collectors.toList());
        List<String> bed = Files.readAllLines(out_bed.toPath()).stream().filter(s -> !s.startsWith("#")).collect(Collectors.toList());;

        assertEquals(vc.stream().filter(v -> v.getAlternateAllele(0).isSingleBreakend()).count(), bed.size());
        // 1 BEDPE per breakend pair
        assertEquals(vc.stream().filter(v -> v.getAlternateAllele(0).isBreakpoint()).count(), 2 * bedpe.size());
    }
    @Test
    @Category(Hg19Tests.class)
    public void should_process_per_sample() throws IOException {
        File out_bed = new File(testFolder.getRoot(), "out.bed");
        File out_bedpe = new File(testFolder.getRoot(), "out.bedpe");
        new GeneratePonBedpe().instanceMain(new String[] {
                "INPUT=" + colo829.getAbsolutePath(),
                "OUTPUT_BED=" + out_bed.getAbsolutePath(),
                "OUTPUT_BEDPE=" + out_bedpe.getAbsolutePath(),
                "NORMAL_ORDINAL=0",
                "R=" + Hg19Tests.findBroadHg19Reference().getAbsolutePath()
        });
        List<String> bedpe = Files.readAllLines(out_bedpe.toPath()).stream().filter(s -> !s.startsWith("#")).collect(Collectors.toList());
        List<String> bed = Files.readAllLines(out_bed.toPath()).stream().filter(s -> !s.startsWith("#")).collect(Collectors.toList());;

        assertEquals(0, bed.size());
        assertEquals(0, bedpe.size());
    }
    @Test
    @Category(Hg19Tests.class)
    public void should_append_to_existing_pon() throws IOException {
        File out_bed = new File(testFolder.getRoot(), "out.bed");
        File out_bedpe = new File(testFolder.getRoot(), "out.bedpe");
        File out2_bed = new File(testFolder.getRoot(), "out2.bed");
        File out2_bedpe = new File(testFolder.getRoot(), "out2.bedpe");
        new GeneratePonBedpe().instanceMain(new String[] {
                "INPUT=" + colo829.getAbsolutePath(),
                "OUTPUT_BED=" + out_bed.getAbsolutePath(),
                "OUTPUT_BEDPE=" + out_bedpe.getAbsolutePath(),
                "NORMAL_ORDINAL=1",
                "R=" + Hg19Tests.findBroadHg19Reference().getAbsolutePath()
        });
        new GeneratePonBedpe().instanceMain(new String[] {
                "INPUT=" + colo829.getAbsolutePath(),
                "OUTPUT_BED=" + out2_bed.getAbsolutePath(),
                "OUTPUT_BEDPE=" + out2_bedpe.getAbsolutePath(),
                "INPUT_BED=" + out_bed.getAbsolutePath(),
                "INPUT_BEDPE=" + out_bedpe.getAbsolutePath(),
                "NORMAL_ORDINAL=1",
                "R=" + Hg19Tests.findBroadHg19Reference().getAbsolutePath()
        });
        List<IdsvVariantContext> vc = getVcf(colo829, getBroadHg19PC(), null);
        List<String> bedpe1 = Files.readAllLines(out_bedpe.toPath()).stream().filter(s -> !s.startsWith("#")).collect(Collectors.toList());
        List<String> bed1 = Files.readAllLines(out_bed.toPath()).stream().filter(s -> !s.startsWith("#")).collect(Collectors.toList());;
        List<String> bedpe2 = Files.readAllLines(out2_bedpe.toPath()).stream().filter(s -> !s.startsWith("#")).collect(Collectors.toList());
        List<String> bed2 = Files.readAllLines(out2_bed.toPath()).stream().filter(s -> !s.startsWith("#")).collect(Collectors.toList());;

        assertEquals(bed1.size(), bed2.size());
        assertEquals(bedpe1.size(), bedpe2.size());

        for (Pair<String, String> p : Streams.zip(bed1.stream(), bed2.stream(), (a, b) -> Pair.create(a, b)).collect(Collectors.toList())) {
            String[] s1 = p.getFirst().split("\t");
            String[] s2 = p.getSecond().split("\t");
            for (int i = 0; i < s1.length; i++) {
                if (i == 4) {
                    assertEquals("1", s1[i]);
                    assertEquals("2", s2[i]);
                } else {
                    assertEquals(s1[i], s2[i]);
                }
            }
        }
        for (Pair<String, String> p : Streams.zip(bedpe1.stream(), bedpe2.stream(), (a, b) -> Pair.create(a, b)).collect(Collectors.toList())) {
            String[] s1 = p.getFirst().split("\t");
            String[] s2 = p.getSecond().split("\t");
            for (int i = 0; i < s1.length; i++) {
                if (i == 7) {
                    assertEquals("1", s1[i]);
                    assertEquals("2", s2[i]);
                } else {
                    assertEquals(s1[i], s2[i]);
                }
            }
        }
    }
    @Test
    @Category(Hg19Tests.class)
    public void should_round_trip() throws IOException {
        File out_bed = new File(testFolder.getRoot(), "out.bed");
        File out_bedpe = new File(testFolder.getRoot(), "out.bedpe");
        File out2_bed = new File(testFolder.getRoot(), "out2.bed");
        File out2_bedpe = new File(testFolder.getRoot(), "out2.bedpe");
        new GeneratePonBedpe().instanceMain(new String[] {
                "INPUT=" + colo829.getAbsolutePath(),
                "OUTPUT_BED=" + out_bed.getAbsolutePath(),
                "OUTPUT_BEDPE=" + out_bedpe.getAbsolutePath(),
                "NORMAL_ORDINAL=1",
                "R=" + Hg19Tests.findBroadHg19Reference().getAbsolutePath()
        });
        List<String> bedpe1 = Files.readAllLines(out_bedpe.toPath()).stream().filter(s -> !s.startsWith("#")).collect(Collectors.toList());
        List<String> bed1 = Files.readAllLines(out_bed.toPath()).stream().filter(s -> !s.startsWith("#")).collect(Collectors.toList());
        new GeneratePonBedpe().instanceMain(new String[] {
                "INPUT=" + colo829.getAbsolutePath(),
                "OUTPUT_BED=" + out2_bed.getAbsolutePath(),
                "OUTPUT_BEDPE=" + out2_bedpe.getAbsolutePath(),
                "INPUT_BED=" + out_bed.getAbsolutePath(),
                "INPUT_BEDPE=" + out_bedpe.getAbsolutePath(),
                "NORMAL_ORDINAL=0",
                "R=" + Hg19Tests.findBroadHg19Reference().getAbsolutePath()
        });
        assertEquals(S(Files.readAllBytes(out_bed.toPath())), S(Files.readAllBytes(out2_bed.toPath())));
        assertEquals(S(Files.readAllBytes(out_bedpe.toPath())), S(Files.readAllBytes(out2_bedpe.toPath())));
    }
}