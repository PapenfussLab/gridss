package gridss;

import au.edu.wehi.idsv.*;
import au.edu.wehi.idsv.configuration.GridssConfiguration;
import au.edu.wehi.idsv.picard.SynchronousReferenceLookupAdapter;
import htsjdk.samtools.SAMFileWriterImpl;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import org.junit.Assert;
import org.junit.Test;
import org.junit.experimental.categories.Category;

import java.io.File;
import java.io.IOException;

import static au.edu.wehi.idsv.TestHelper.getConfig;

public class SanityCheckEvidenceTest {
    //@Test // fixed in preprocessing
    @Category(EColiTests.class)
    public void issue278_inconsistent_split_read() throws IOException {
        File ref = ReferenceTests.findReference("Escherichia_coli_bl21_de3_.ASM956v1.dna.toplevel.fa");
        SynchronousReferenceLookupAdapter reflookup = new SynchronousReferenceLookupAdapter(new IndexedFastaSequenceFile(ref));
        File input = new File("src/test/resources/sanity_failure_debug/bl21_de3_.ASM956v1_S2.bam");
        ProcessingContext pc = new ProcessingContext(new FileSystemContext(input.getParentFile(), input.getParentFile(), SAMFileWriterImpl.getDefaultMaxRecordsInRam()), ref, reflookup, null, getConfig());
        SAMEvidenceSource source = new SAMEvidenceSource(pc, input, null, 0);
        SanityCheckEvidence sce = new SanityCheckEvidence();
        Assert.assertEquals(0, sce.sanityCheck(source));
    }
    @Test
    @Category(EColiTests.class)
    public void issue278_inconsistent_assembly_scoring() throws IOException {
        File ref = ReferenceTests.findReference("Escherichia_coli_bl21_de3_.ASM956v1.dna.toplevel.fa");
        SynchronousReferenceLookupAdapter reflookup = new SynchronousReferenceLookupAdapter(new IndexedFastaSequenceFile(ref));
        File input = new File("src/test/resources/sanity_failure_debug/bl21_de3_.ASM956v1_assembly.bam");
        GridssConfiguration config = getConfig();
        config.hashEvidenceID = false;
        ProcessingContext pc = new ProcessingContext(new FileSystemContext(input.getParentFile(), input.getParentFile(), SAMFileWriterImpl.getDefaultMaxRecordsInRam()), ref, reflookup, null, config);
        SAMEvidenceSource source = new SAMEvidenceSource(pc, input, null, 0);
        SanityCheckEvidence sce = new SanityCheckEvidence();
        Assert.assertEquals(0, sce.sanityCheck(source));
    }
    @Test
    @Category(Hg19Tests.class)
    public void colo829_1_inconsistent_split_read_scoring() throws IOException {
        File ref = ReferenceTests.findReference("Homo_sapiens.GRCh37.GATK.illumina.fasta");
        SynchronousReferenceLookupAdapter reflookup = new SynchronousReferenceLookupAdapter(new IndexedFastaSequenceFile(ref));
        File input = new File("src/test/resources/sanity_failure_debug/colo829_1/COLO829v001R_dedup.realigned.bam");
        GridssConfiguration config = getConfig();
        config.hashEvidenceID = false;
        ProcessingContext pc = new ProcessingContext(new FileSystemContext(input.getParentFile(), input.getParentFile(), SAMFileWriterImpl.getDefaultMaxRecordsInRam()), ref, reflookup, null, config);
        SAMEvidenceSource source = new SAMEvidenceSource(pc, input, null, 0);
        SanityCheckEvidence sce = new SanityCheckEvidence();
        Assert.assertEquals(0, sce.sanityCheck(source));
    }
}