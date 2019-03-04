package gridss;

import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Random;

@CommandLineProgramProperties(
        summary = "Simple utility program to test robustness of deletion and duplication calls.",
        oneLineSummary = "Simple utility program to test robustness of deletion and duplication calls",
        programGroup = gridss.cmdline.programgroups.Benchmarking.class
)
public class SyntheticDelDupSimulator extends CommandLineProgram {
    private static final Log log = Log.getInstance(AnnotateInexactHomologyBedpe.class);
    @Argument(doc="Reference genome")
    public File REFERENCE_OUTPUT;
    @Argument(doc="Sequence of genome with variants inserted")
    public File VARIANT_OUTPUT;
    @Argument(doc="Bases around the events")
    public int FLANKING_BASES = 1000;
    @Argument(doc="Event sizes. Negative indicates deletion, positive indicates duplication.")
    public List<Integer> SIZE;
    @Argument(doc="Event homologies. Positve indicates homology around the site, negative indicates additional inserted sequence")
    public List<Integer> INSHOM;

    public static void main(String[] argv) {
        System.exit(new SyntheticDelDupSimulator().instanceMain(argv));
    }

    @Override
    protected int doWork() {
        Random rng = new Random(0);
        try (BufferedWriter wref = new BufferedWriter(new FileWriter(REFERENCE_OUTPUT))) {
            try (BufferedWriter wvar = new BufferedWriter(new FileWriter(VARIANT_OUTPUT))) {
                for (int i = 0; i < SIZE.size(); i++) {
                    int len = SIZE.get(i);
                    int inslen = INSHOM.get(i);
                    int hombases = Math.max(0, inslen);
                    int insbases = -Math.min(0, inslen);
                    String leftFlank = getRandomBases(rng, FLANKING_BASES);
                    String rightFlank = getRandomBases(rng, FLANKING_BASES);
                    String ins = getRandomBases(rng, insbases);
                    String hom = getRandomBases(rng, hombases);
                    String refSeq;
                    String varSeq;
                    String type;
                    if (len < 0) {
                        type = "DEL";
                        len *= -1;
                        String deletedBases = getRandomBases(rng, len - hombases);
                        // deletion could be (hom + deletedBased), or (deletedBases + hom)
                        refSeq = leftFlank + hom + deletedBases + hom + rightFlank;
                        varSeq = leftFlank + hom + ins + rightFlank;
                    } else {
                        type = "DUP";
                        // DUP
                        String dupBases = getRandomBases(rng, len);
                        refSeq = leftFlank + hom + dupBases + hom + rightFlank;
                        varSeq = leftFlank + hom + dupBases + hom + ins + dupBases + hom + rightFlank;
                    }
                    String header = ">" + getContigName(i, FLANKING_BASES, hombases, insbases, type, len);
                    wref.write(header + "\n" + refSeq + "\n");
                    wvar.write(header + "\n" + varSeq + "\n");
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return 0;
    }
    private static char[] bases = new char[] {'A','C', 'G', 'T' };
    private static String getRandomBases(Random rng, int len) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < len; i++) {
            sb.append(bases[rng.nextInt(bases.length)]);
        }
        return sb.toString();
    }
    private static String getContigName(int ordinal, int flanking, int homlen, int inslen, String type, int size) {
        return String.format("contig%d_%s%d_hom%d_ins%d", ordinal, type, size, homlen, inslen);
    }
}
