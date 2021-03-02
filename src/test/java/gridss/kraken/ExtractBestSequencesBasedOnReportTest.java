package gridss.kraken;

import au.edu.wehi.idsv.IntermediateFilesTest;
import au.edu.wehi.idsv.ncbi.TaxonomyLevel;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

// for f in $(find src/test/resources/kraken        -name '*.fa') ; do grep ">" $f | cut -b 2- | awk '{ split($0,a,"|"); split($0,b," ") ; print b[1] "\t" a[2]  } ' >> $(dirname $f)/seqid2taxid.map ; done
// for f in $(find src/test/resources/virusbreakend -name '*.fa') ; do grep ">" $f | cut -b 2- | awk '{ split($0,a,"|"); split($0,b," ") ; print b[1] "\t" a[2]  } ' >> $(dirname $f)/seqid2taxid.map ; done
public class ExtractBestSequencesBasedOnReportTest extends IntermediateFilesTest {
    @Test
    public void regression_should_extract_from_multiple_fasta() throws IOException {
        ExtractBestSequencesBasedOnReport cmd = new ExtractBestSequencesBasedOnReport();
        cmd.INPUT_KRAKEN2_REPORT = new File("src/test/resources/kraken/multiple_dictionaries/all_subspecies.txt");
        cmd.KRAKEN_REFERENCES = ImmutableList.of(
                new File("src/test/resources/kraken/multiple_dictionaries/kraken.fa"),
                new File("src/test/resources/kraken/multiple_dictionaries/virushostdb.fa")
        );
        cmd.NCBI_NODES_DMP = new File("src/test/resources/kraken/multiple_dictionaries/matching.nodes.dmp");
        cmd.SEQID2TAXID_MAP = new File("src/test/resources/kraken/multiple_dictionaries/seqid2taxid.map");
        cmd.OUTPUT = new File(output.toString() + ".fa");;
        cmd.REPORT_OUTPUT = new File(output.toString() + ".report.txt");
        cmd.SUMMARY_REPORT_OUTPUT = new File(output.toString() + ".summary.report.txt");
        cmd.TAXONOMY_IDS = ImmutableList.of(10239);
        cmd.doWork();

        Assert.assertTrue(cmd.OUTPUT.exists());
        Assert.assertTrue(cmd.REPORT_OUTPUT.exists());
        Assert.assertTrue(cmd.SUMMARY_REPORT_OUTPUT.exists());

        List<String> out = Files.readAllLines(cmd.OUTPUT.toPath());
        List<String> report = Files.readAllLines(cmd.REPORT_OUTPUT.toPath());
        List<String> summary = Files.readAllLines(cmd.SUMMARY_REPORT_OUTPUT.toPath());
        Assert.assertEquals(1, summary.size());
        Assert.assertTrue(summary.get(0).contains("papillomavirus type 45"));
        Assert.assertTrue(out.get(0).startsWith(">kraken:taxid|10593|X74479"));

    }
    @Test
    public void should_not_extract_all_strains() throws IOException {
        ExtractBestSequencesBasedOnReport cmd = new ExtractBestSequencesBasedOnReport();
        cmd.INPUT_KRAKEN2_REPORT = new File("src/test/resources/kraken/multistrain/26T.virusbreakend.vcf.kraken2.report.all.txt");
        cmd.KRAKEN_REFERENCES = ImmutableList.of(
                new File("src/test/resources/kraken/multistrain/26T.virusbreakend.vcf.viral.fa")
        );
        //grep -f <(cut -f 5 < 26T.virusbreakend.vcf.kraken2.report.all.txt | sed 's/^/^/' | sed 's/$/    /') nodes.dmp > relevant_nodes.dmp
        cmd.NCBI_NODES_DMP = new File("src/test/resources/kraken/multistrain/relevant_nodes.dmp");
        cmd.SEQID2TAXID_MAP = new File("src/test/resources/kraken/multistrain/seqid2taxid.map");
        cmd.OUTPUT = new File(output.toString() + ".fa");;
        cmd.REPORT_OUTPUT = new File(output.toString() + ".report.txt");
        cmd.SUMMARY_REPORT_OUTPUT = new File(output.toString() + ".summary.report.txt");
        cmd.TAXONOMY_IDS = ImmutableList.of(10239);
        cmd.MIN_SUPPORTING_READS = 1;
        cmd.doWork();

        Assert.assertTrue(cmd.OUTPUT.exists());
        ReferenceSequenceFile ref = new FastaSequenceFile(cmd.OUTPUT, true);
        List<ReferenceSequence> refList = new ArrayList<>();
        ReferenceSequence rs = ref.nextSequence();
        while (rs != null) {
            refList.add(rs);
            rs = ref.nextSequence();
        }
        Assert.assertEquals(1, refList.size());
    }
    @Test
    public void should_output_one_reference_per_genus() throws IOException {
        // Script to extract nodes.dmp and reference sequences of interest:
        // grep ~/dev/virusbreakend/virusbreakenddb/taxonomy/nodes.dmp -f <(cut -f 5 *.report.all.txt | sed 's/$/\t/' | sed 's/^/^/') > relevant_nodes.dmp
        // samtools faidx ~/dev/virusbreakend/virusbreakenddb/library/added/*.fna $(grep ~/dev/virusbreakend/virusbreakenddb/library/added/*.fna -f <(cut -f 5 *.report.all.txt | sed 's/$/[|]/' | sed 's/^/>kraken:taxid[|]/') | cut -f 1 -d ' ' | cut -b 2- | sed 's/$/:1-100000/' | tr '\n' " ") | sed 's/:1-100000//' > relevant_ref.fa
        // samtools faidx relevant_ref.fa && picard CreateSequenceDictionary R=relevant_ref.fa
        ExtractBestSequencesBasedOnReport cmd = new ExtractBestSequencesBasedOnReport();
        cmd.INPUT_KRAKEN2_REPORT = new File("src/test/resources/kraken/genus/hbv.kraken2.report.all.txt");
        cmd.KRAKEN_REFERENCES = ImmutableList.of(
                new File("src/test/resources/kraken/genus/relevant_ref.fa")
        );
        //grep -f <(cut -f 5 < 26T.virusbreakend.vcf.kraken2.report.all.txt | sed 's/^/^/' | sed 's/$/    /') nodes.dmp > relevant_nodes.dmp
        cmd.NCBI_NODES_DMP = new File("src/test/resources/kraken/genus/relevant_nodes.dmp");
        cmd.SEQID2TAXID_MAP = new File("src/test/resources/kraken/genus/seqid2taxid.map");
        cmd.OUTPUT = new File(output.toString() + ".fa");;
        cmd.REPORT_OUTPUT = new File(output.toString() + ".report.txt");
        cmd.SUMMARY_REPORT_OUTPUT = new File(output.toString() + ".summary.report.txt");
        cmd.TAXONOMY_IDS = ImmutableList.of(10239);
        cmd.MIN_SUPPORTING_READS = 1;
        cmd.doWork();

        Assert.assertTrue(cmd.OUTPUT.exists());
        ReferenceSequenceFile ref = new FastaSequenceFile(cmd.OUTPUT, true);
        List<ReferenceSequence> refList = new ArrayList<>();
        ReferenceSequence rs = ref.nextSequence();
        while (rs != null) {
            refList.add(rs);
            rs = ref.nextSequence();
        }
        Assert.assertEquals(1, refList.size());
    }
    @Test
    public void should_output_per_species() throws IOException {
        // Script to extract nodes.dmp and reference sequences of interest:
        // grep ~/dev/virusbreakend/virusbreakenddb/taxonomy/nodes.dmp -f <(cut -f 5 *.report.all.txt | sed 's/$/\t/' | sed 's/^/^/') > relevant_nodes.dmp
        // samtools faidx ~/dev/virusbreakend/virusbreakenddb/library/added/*.fna $(grep ~/dev/virusbreakend/virusbreakenddb/library/added/*.fna -f <(cut -f 5 *.report.all.txt | sed 's/$/[|]/' | sed 's/^/>kraken:taxid[|]/') | cut -f 1 -d ' ' | cut -b 2- | sed 's/$/:1-100000/' | tr '\n' " ") | sed 's/:1-100000//' > relevant_ref.fa
        // samtools faidx relevant_ref.fa && picard CreateSequenceDictionary R=relevant_ref.fa
        ExtractBestSequencesBasedOnReport cmd = new ExtractBestSequencesBasedOnReport();
        String datadir = "src/test/resources/kraken/genus/";
        cmd.INPUT_KRAKEN2_REPORT = new File(datadir + "hbv.kraken2.report.all.txt");
        cmd.NCBI_NODES_DMP = new File(datadir + "relevant_nodes.dmp");
        cmd.SEQID2TAXID_MAP = new File(datadir + "seqid2taxid.map");
        cmd.KRAKEN_REFERENCES = ImmutableList.of(
                new File(datadir + "relevant_ref.fa")
        );
        cmd.OUTPUT = new File(output.toString() + ".fa");;
        cmd.REPORT_OUTPUT = new File(output.toString() + ".report.txt");
        cmd.SUMMARY_REPORT_OUTPUT = new File(output.toString() + ".summary.report.txt");
        cmd.TAXONOMY_IDS = ImmutableList.of(10239);
        cmd.MIN_SUPPORTING_READS = 1;
        cmd.TAXONOMIC_DEDUPLICATION_LEVEL = null;
        cmd.doWork();

        Assert.assertTrue(cmd.OUTPUT.exists());
        ReferenceSequenceFile ref = new FastaSequenceFile(cmd.OUTPUT, true);
        List<ReferenceSequence> refList = new ArrayList<>();
        ReferenceSequence rs = ref.nextSequence();
        while (rs != null) {
            refList.add(rs);
            rs = ref.nextSequence();
        }
        Assert.assertEquals(15, refList.size());
    }
    @Test
    public void herpes_should_not_mask_HBV() throws IOException {
        ExtractBestSequencesBasedOnReport cmd = new ExtractBestSequencesBasedOnReport();
        String datadir = "src/test/resources/kraken/genus/herpes_should_not_mask_HBV/";
        cmd.INPUT_KRAKEN2_REPORT = new File(datadir + "herpes_HBV.kraken2.report.all.txt");
        cmd.NCBI_NODES_DMP = new File(datadir + "relevant_nodes.dmp");
        cmd.SEQID2TAXID_MAP = new File(datadir + "seqid2taxid.map");
        cmd.KRAKEN_REFERENCES = ImmutableList.of(
                new File(datadir + "relevant_ref.fa")
        );
        cmd.OUTPUT = new File(output.toString() + ".fa");;
        cmd.REPORT_OUTPUT = new File(output.toString() + ".report.txt");
        cmd.SUMMARY_REPORT_OUTPUT = new File(output.toString() + ".summary.report.txt");
        cmd.MIN_SUPPORTING_READS = 50;
        cmd.doWork();

        Assert.assertTrue(cmd.OUTPUT.exists());
        ReferenceSequenceFile ref = new FastaSequenceFile(cmd.OUTPUT, true);
        List<ReferenceSequence> refList = new ArrayList<>();
        ReferenceSequence rs = ref.nextSequence();
        while (rs != null) {
            refList.add(rs);
            rs = ref.nextSequence();
        }
        Assert.assertEquals(2, refList.size());
    }
    @Test
    public void summary_output_should_traverse_species_genus() throws IOException {
        ExtractBestSequencesBasedOnReport cmd = new ExtractBestSequencesBasedOnReport();
        cmd.INPUT_KRAKEN2_REPORT = new File("src/test/resources/kraken/multistrain/26T.virusbreakend.vcf.kraken2.report.all.txt");
        cmd.KRAKEN_REFERENCES = ImmutableList.of(
                new File("src/test/resources/kraken/multistrain/26T.virusbreakend.vcf.viral.fa")
        );
        cmd.NCBI_NODES_DMP = new File("src/test/resources/kraken/multistrain/relevant_nodes.dmp");
        cmd.SEQID2TAXID_MAP = new File("src/test/resources/kraken/multistrain/seqid2taxid.map");
        cmd.OUTPUT = new File(output.toString() + ".fa");;
        cmd.REPORT_OUTPUT = new File(output.toString() + ".report.txt");
        cmd.SUMMARY_REPORT_OUTPUT = new File(output.toString() + ".summary.report.txt");
        cmd.SUMMARY_OUTPUT = new File(output.toString() + ".summary.tsv");
        cmd.TAXONOMY_IDS = ImmutableList.of(10239);
        cmd.MIN_SUPPORTING_READS = 1;
        cmd.doWork();

        List<String> summary = Files.readAllLines(cmd.SUMMARY_OUTPUT.toPath());
        Assert.assertEquals("10405\tOrthohepadnavirus\t858\t10407\tHepatitis B virus\t856\t489466\tHBV genotype C\t277\t276\tkraken:taxid|489466|AB540583", summary.get(1));
    }
    //@Test
    @Ignore // Underlying issue was that HPV-45 isn't in ReqSeq and the first half of VIRUSBreakend was run against a larger database
    public void regression_should_call_hpv45() throws IOException {
        // Script to extract nodes.dmp and reference sequences of interest:
        // grep ~/dev/virusbreakend/virusbreakenddb/taxonomy/nodes.dmp -f <(cut -f 5 *.report.all.txt | sed 's/$/\t/' | sed 's/^/^/') > relevant_nodes.dmp
        // samtools faidx ~/dev/virusbreakend/virusbreakenddb/library/viral/*.fna $(grep ~/dev/virusbreakend/virusbreakenddb/library/viral/*.fna -f <(cut -f 5 *.report.all.txt | sed 's/$/[|]/' | sed 's/^/>kraken:taxid[|]/') | cut -f 1 -d ' ' | cut -b 2- | sed 's/$/:1-1000/' | tr '\n' " ") | sed 's/:1-1000//' > relevant_ref.fa
        // samtools faidx relevant_ref.fa && picard CreateSequenceDictionary R=relevant_ref.fa
        ExtractBestSequencesBasedOnReport cmd = new ExtractBestSequencesBasedOnReport();
        String datadir = "src/test/resources/virusbreakend/hpv45/";
        cmd.INPUT_KRAKEN2_REPORT = new File(datadir + "hpv45.virusbreakend.vcf.kraken2.report.all.txt");
        cmd.NCBI_NODES_DMP = new File(datadir + "relevant_nodes.dmp");
        cmd.SEQID2TAXID_MAP = new File(datadir + "seqid2taxid.map");
        cmd.KRAKEN_REFERENCES = ImmutableList.of(new File(datadir + "relevant_ref.fa"));
        cmd.OUTPUT = new File(output.toString() + ".fa");;
        cmd.REPORT_OUTPUT = new File(output.toString() + ".report.txt");
        cmd.SUMMARY_REPORT_OUTPUT = new File(output.toString() + ".summary.report.txt");
        cmd.TAXONOMY_IDS = ImmutableList.of(10239);
        cmd.TAXONOMIC_DEDUPLICATION_LEVEL = TaxonomyLevel.Genus;
        cmd.MIN_SUPPORTING_READS = 50;
        cmd.doWork();
        List<String> summary = Files.readAllLines(cmd.SUMMARY_REPORT_OUTPUT.toPath());
        Assert.assertTrue(summary.get(0).contains("Human papillomavirus type 45"));
    }
}