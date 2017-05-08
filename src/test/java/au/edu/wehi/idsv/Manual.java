package au.edu.wehi.idsv;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.commons.configuration.ConfigurationException;
import org.junit.Assert;
import org.junit.Test;
import org.junit.experimental.categories.Category;

import au.edu.wehi.idsv.configuration.GridssConfiguration;
import au.edu.wehi.idsv.picard.SynchronousReferenceLookupAdapter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

/**
 * Ad-hoc debugging tests
 * @author Daniel Cameron
 *
 */
public class Manual extends IntermediateFilesTest {
	/**
	 * Test our iterators are behaving correctly
	 * @throws IOException 
	 */
	//@Test
	@Category(Hg19Tests.class)
	public void debug778sorting() throws ConfigurationException, IOException {
		ProcessingContext pc = new ProcessingContext(
			new FileSystemContext(new File("W:\\778\\idsv"), new File("W:\\778\\idsv"), 500000), Hg19Tests.findHg19Reference(), null, new ArrayList<Header>(),
			new GridssConfiguration());
		List<SAMEvidenceSource> samEvidence = new ArrayList<SAMEvidenceSource>();
		for (String s : new String[] {
				"W:\\778\\DNA_778_HiSeq_35nt_PE1_bt2_s_rg_cleaned.bam",
				"W:\\778\\DNA_778_IL_35nt_PE1_bt2_s_rg_cleaned.bam",
				"W:\\778\\DNA_778_IL_75nt_PE1_bt2_s_rg_cleaned.bam",
				"W:\\778\\DNA_778_PM_75nt_PE1_bt2_s_rg_cleaned.bam",
				"W:\\778\\DNA_778_PM_lane1_100nt_PE1_bt2_s_rg_cleaned.bam",
				"W:\\778\\DNA_778_PM_lane2_100nt_PE1_bt2_s_rg_cleaned.bam",
			}) {
			SAMEvidenceSource ses = new SAMEvidenceSource(pc, new File(s), null, 0);
			//Iterator<DirectedEvidence> it = ses.iterator(true, true, true);
			//while (it.hasNext()) {
			//	it.next();
			//}
			samEvidence.add(ses);
		}
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, samEvidence, new File("W:\778\\idsv\\778.vcf.idsv.working"));
		aes.assembleBreakends(null);
		//Iterator<SAMRecordAssemblyEvidence> it = aes.iterator(true, true);
		//while (it.hasNext()) {
		//	it.next();
		//}
		Iterator<DirectedEvidence> allIt = SAMEvidenceSource.mergedIterator(samEvidence, false);
		while (allIt.hasNext()) {
			allIt.next();
		}
	}
	//@Test
	@Category(Hg19Tests.class)
	public void assembly_scoring_should_be_symmetrical() throws FileNotFoundException {
		File sam = new File("W:/debug/test.bam");
		File ref = Hg19Tests.findHg19Reference();
		ProcessingContext pc = new ProcessingContext(getFSContext(), ref, new SynchronousReferenceLookupAdapter(new IndexedFastaSequenceFile(ref)), null, getConfig());
		SAMEvidenceSource ses = new SAMEvidenceSource(pc, sam, null, 0);
		List<SAMRecord> reads = getRecords(new File("W:/debug/test.bam.gridss.working/test.bam.sv.bam"));
		List<SplitReadEvidence> list = ses.iterator().stream()
			.filter(e -> e instanceof SplitReadEvidence)
			.map(e -> (SplitReadEvidence)e)
			.filter(e -> e.getSAMRecord().getReadName().equals("variant.chr12-art1000553"))
			.collect(Collectors.toList());
		for (SplitReadEvidence e : list) {
			e.getBreakpointQual();
		}
	}
	@Test
	@Category(Hg38Tests.class)
	public void should_not_include_concordant_read_pairs() throws FileNotFoundException {
		File sam = new File("W:/hg38debug/D1_S1_L001_R1_001.bam.gridss.working/D1_S1_L001_R1_001.bam");
		File ref = Hg38Tests.findHg38Reference();
		ProcessingContext pc = new ProcessingContext(getFSContext(), ref, new SynchronousReferenceLookupAdapter(new IndexedFastaSequenceFile(ref)), null, getConfig());
		SAMEvidenceSource ses = new SAMEvidenceSource(pc, sam, null, 0, 0.995);
		List<DirectedEvidence> list = ses.iterator().stream()
			.collect(Collectors.toList());
		List<DirectedEvidence> reads = list.stream().filter(e -> e instanceof DiscordantReadPair)
			.map(e -> (DiscordantReadPair)e)
			.filter(e -> e.getLocalledMappedRead().getReadName().equals("D00360:95:H2YWMBCXX:1:2107:5758:68093"))
			.collect(Collectors.toList());
		Assert.assertEquals(0, reads.size());
	}
}
