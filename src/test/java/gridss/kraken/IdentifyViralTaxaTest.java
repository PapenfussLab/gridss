package gridss.kraken;

import au.edu.wehi.VirusbreakenddbTests;
import au.edu.wehi.idsv.IntermediateFilesTest;
import au.edu.wehi.idsv.ncbi.TaxonomyLevel;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import org.apache.commons.lang3.tuple.Pair;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

// for f in $(find src/test/resources/kraken        -name '*.fa') ; do grep ">" $f | cut -b 2- | awk '{ split($0,a,"|"); split($0,b," ") ; print b[1] "\t" a[2]  } ' >> $(dirname $f)/seqid2taxid.map ; done
// for f in $(find src/test/resources/virusbreakend -name '*.fa') ; do grep ">" $f | cut -b 2- | awk '{ split($0,a,"|"); split($0,b," ") ; print b[1] "\t" a[2]  } ' >> $(dirname $f)/seqid2taxid.map ; done
public class IdentifyViralTaxaTest extends IntermediateFilesTest {
    private File report;
    private File summaryReport;
    private File intermediateTsv;
    @Before
    @Override
    public void setup() throws IOException {
        super.setup();
        report = new File(output.toString() + ".report.txt");
        summaryReport = new File(output.toString() + ".summary.report.txt");
        intermediateTsv = new File(output.toString() + ".summary.tsv");
    }
    private IdentifyViralTaxa setup(
            File krakenReport,
            List<File> ref,
            File nodesdmp,
            File seqidmap,
            List<Integer> taxids,
            int minSupportingReads,
            TaxonomyLevel level) {
        IdentifyViralTaxa cmd1 = new IdentifyViralTaxa();
        cmd1.INPUT_KRAKEN2_REPORT = krakenReport;
        cmd1.KRAKEN_REFERENCES = ref;
        cmd1.NCBI_NODES_DMP = nodesdmp;
        cmd1.SEQID2TAXID_MAP = seqidmap;
        cmd1.OUTPUT = intermediateTsv;
        cmd1.SUMMARY_REPORT_OUTPUT = summaryReport;
        cmd1.REPORT_OUTPUT = report;
        cmd1.MIN_SUPPORTING_READS = minSupportingReads;
        cmd1.TAXONOMIC_DEDUPLICATION_LEVEL = level;
        if (taxids != null) cmd1.TAXONOMY_IDS = taxids;
        cmd1.doWork();
        return cmd1;
    }
    @Test
    public void regression_should_extract_from_multiple_fasta() throws IOException {
        String datadir = "src/test/resources/kraken/multiple_dictionaries/";
        IdentifyViralTaxa cmd = setup(
                new File(datadir + "all_subspecies.txt"),
                ImmutableList.of(
                        new File(datadir + "kraken.fa"),
                        new File(datadir + "virushostdb.fa")
                ),
                new File(datadir + "matching.nodes.dmp"),
                new File(datadir + "seqid2taxid.map"),
                ImmutableList.of(10239),
                50,
                TaxonomyLevel.Genus);
        Assert.assertTrue(cmd.OUTPUT.exists());
        Assert.assertTrue(cmd.REPORT_OUTPUT.exists());
        Assert.assertTrue(cmd.SUMMARY_REPORT_OUTPUT.exists());

        List<String> out = Files.readAllLines(cmd.OUTPUT.toPath());
        List<String> report = Files.readAllLines(cmd.REPORT_OUTPUT.toPath());
        List<String> summary = Files.readAllLines(cmd.SUMMARY_REPORT_OUTPUT.toPath());
        Assert.assertTrue(summary.get(0).contains("papillomavirus type 45"));
        Assert.assertTrue(report.size() > summary.size());
        Assert.assertTrue(out.get(1).contains("papillomavirus type 45"));
        Assert.assertTrue(out.get(1).contains("10593"));
    }
    @Test
    public void should_not_extract_all_strains() throws IOException {
        String datadir = "src/test/resources/kraken/multistrain/";
        IdentifyViralTaxa cmd = setup(
                new File(datadir + "26T.virusbreakend.vcf.kraken2.report.all.txt"),
                ImmutableList.of(new File(datadir + "26T.virusbreakend.vcf.viral.fa")),
                new File(datadir + "relevant_nodes.dmp"),
                new File(datadir + "seqid2taxid.map"),
                ImmutableList.of(10239),
                1,
                TaxonomyLevel.Genus);

        List<String> out = Files.readAllLines(cmd.OUTPUT.toPath());
        Assert.assertEquals(1+1, out.size()); // 1 + header
    }
    @Test
    public void should_output_one_reference_per_genus() throws IOException {
        // Script to extract nodes.dmp and reference sequences of interest:
        // grep ~/dev/virusbreakend/virusbreakenddb/taxonomy/nodes.dmp -f <(cut -f 5 *.report.all.txt | sed 's/$/\t/' | sed 's/^/^/') > relevant_nodes.dmp
        // samtools faidx ~/dev/virusbreakend/virusbreakenddb/library/added/*.fna $(grep ~/dev/virusbreakend/virusbreakenddb/library/added/*.fna -f <(cut -f 5 *.report.all.txt | sed 's/$/[|]/' | sed 's/^/>kraken:taxid[|]/') | cut -f 1 -d ' ' | cut -b 2- | sed 's/$/:1-100000/' | tr '\n' " ") | sed 's/:1-100000//' > relevant_ref.fa
        // samtools faidx relevant_ref.fa && picard CreateSequenceDictionary R=relevant_ref.fa
        String datadir = "src/test/resources/kraken/genus/";
        IdentifyViralTaxa cmd = setup(
            new File(datadir + "hbv.kraken2.report.all.txt"),
            ImmutableList.of(new File(datadir + "relevant_ref.fa")),
            new File(datadir + "relevant_nodes.dmp"),
            new File(datadir + "seqid2taxid.map"),
            ImmutableList.of(10239),
            50,
            TaxonomyLevel.Genus);
        //grep -f <(cut -f 5 < 26T.virusbreakend.vcf.kraken2.report.all.txt | sed 's/^/^/' | sed 's/$/    /') nodes.dmp > relevant_nodes.dmp

        List<String> output = Files.readAllLines(cmd.OUTPUT.toPath());
        Assert.assertEquals(1+1, output.size()); // 1 + header
    }
    @Test
    public void should_output_per_species() throws IOException {
        // Script to extract nodes.dmp and reference sequences of interest:
        // grep ~/dev/virusbreakend/virusbreakenddb/taxonomy/nodes.dmp -f <(cut -f 5 *.report.all.txt | sed 's/$/\t/' | sed 's/^/^/') > relevant_nodes.dmp
        // samtools faidx ~/dev/virusbreakend/virusbreakenddb/library/added/*.fna $(grep ~/dev/virusbreakend/virusbreakenddb/library/added/*.fna -f <(cut -f 5 *.report.all.txt | sed 's/$/[|]/' | sed 's/^/>kraken:taxid[|]/') | cut -f 1 -d ' ' | cut -b 2- | sed 's/$/:1-100000/' | tr '\n' " ") | sed 's/:1-100000//' > relevant_ref.fa
        // samtools faidx relevant_ref.fa && picard CreateSequenceDictionary R=relevant_ref.fa
        String datadir = "src/test/resources/kraken/genus/";
        IdentifyViralTaxa cmd = setup(
            new File(datadir + "hbv.kraken2.report.all.txt"),
            ImmutableList.of(new File(datadir + "relevant_ref.fa")),
            new File(datadir + "relevant_nodes.dmp"),
            new File(datadir + "seqid2taxid.map"),
            ImmutableList.of(10239),
            1,
            null);
        List<String> output = Files.readAllLines(cmd.OUTPUT.toPath());
        Assert.assertEquals(15+1, output.size());
    }
    @Test
    public void herpes_should_not_mask_HBV() throws IOException {
        String datadir = "src/test/resources/kraken/genus/herpes_should_not_mask_HBV/";
        IdentifyViralTaxa cmd = setup(
                new File(datadir + "herpes_HBV.kraken2.report.all.txt"),
                ImmutableList.of(new File(datadir + "relevant_ref.fa")),
                new File(datadir + "relevant_nodes.dmp"),
                new File(datadir + "seqid2taxid.map"),
                ImmutableList.of(10239),
                50,
                TaxonomyLevel.Genus);
        List<String> output = Files.readAllLines(cmd.OUTPUT.toPath());
        Assert.assertEquals(2+1, output.size());
    }
    @Test
    public void summary_output_should_traverse_species_genus() throws IOException {
        IdentifyViralTaxa cmd = new IdentifyViralTaxa();
        cmd.INPUT_KRAKEN2_REPORT = new File("src/test/resources/kraken/multistrain/26T.virusbreakend.vcf.kraken2.report.all.txt");
        cmd.KRAKEN_REFERENCES = ImmutableList.of(
                new File("src/test/resources/kraken/multistrain/26T.virusbreakend.vcf.viral.fa")
        );
        cmd.NCBI_NODES_DMP = new File("src/test/resources/kraken/multistrain/relevant_nodes.dmp");
        cmd.SEQID2TAXID_MAP = new File("src/test/resources/kraken/multistrain/seqid2taxid.map");
        cmd.OUTPUT = new File(output.toString() + ".fa");;
        cmd.REPORT_OUTPUT = new File(output.toString() + ".report.txt");
        cmd.SUMMARY_REPORT_OUTPUT = new File(output.toString() + ".summary.report.txt");
        cmd.TAXONOMY_IDS = ImmutableList.of(10239);
        cmd.MIN_SUPPORTING_READS = 1;
        cmd.doWork();

        List<String> summary = Files.readAllLines(cmd.OUTPUT.toPath());
        Assert.assertEquals("10405\tOrthohepadnavirus\t858\t10407\tHepatitis B virus\t856\t489466\tHBV genotype C\t277\t276", summary.get(1));
    }
    /* // root cause was blank lines in the host filter in the virusbreakend.sh script
    @Test
    @Category(VirusbreakenddbTests.class)
    public void should_filter_to_human() throws IOException {
        IdentifyViralTaxa cmd = new IdentifyViralTaxa();
        cmd.INPUT_KRAKEN2_REPORT = new File("src/test/resources/kraken/");
        cmd.KRAKEN_REFERENCES = ImmutableList.of(VirusbreakenddbTests.findVirusbreakendFile("library.fna"));
        cmd.NCBI_NODES_DMP = VirusbreakenddbTests.getFullNodesDmp();
        cmd.SEQID2TAXID_MAP = VirusbreakenddbTests.getSeqid2taxidMap();
        cmd.OUTPUT = new File(output.toString() + ".fa");;
        cmd.REPORT_OUTPUT = new File(output.toString() + ".report.txt");
        cmd.SUMMARY_REPORT_OUTPUT = new File(output.toString() + ".summary.report.txt");
        cmd.TAXONOMY_IDS = ImmutableList.of(10239);
        cmd.MIN_SUPPORTING_READS = 1;
        cmd.doWork();
    }
    */
}