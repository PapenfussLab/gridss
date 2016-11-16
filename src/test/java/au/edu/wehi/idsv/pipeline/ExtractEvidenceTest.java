package au.edu.wehi.idsv.pipeline;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.EnumSet;
import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.IntermediateFilesTest;
import au.edu.wehi.idsv.ProcessStep;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMEvidenceSource;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;


public class ExtractEvidenceTest extends IntermediateFilesTest {
	@Test
	public void should_generate_rp_sc_mate_fq() {
		createInput();
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, 0);
		ExtractEvidence e = new ExtractEvidence(getCommandlineContext(), source); e.process(EnumSet.allOf(ProcessStep.class)); e.close();
		assertEquals(0, getRP(source).size());
		assertEquals(0, getSC(source).size());
		assertEquals(0, getMate(source).size());
		assertEquals(0, getFastqRecords(source).size());
	}
	@Test
	public void should_generate_rp_sc_mate_fq_per_chr() {
		ProcessingContext context = getCommandlineContext(true);
		createInput();
		SAMEvidenceSource source = new SAMEvidenceSource(context, input, 0);
		ExtractEvidence e = new ExtractEvidence(context, source); e.process(EnumSet.allOf(ProcessStep.class)); e.close();
		assertEquals(0, new PerChr(context).getRP(source).size());
		assertEquals(0, new PerChr(context).getSC(source).size());
		assertEquals(0, new PerChr(context).getMate(source).size());
		assertEquals(0, new PerChr(context).getFastqRecords(source).size());
	}
	@Test
	public void sc_should_be_located_in_sc_bam_for_chr() {
		ProcessingContext context = getCommandlineContext(true);
		createInput(Read(1, 1, "50M50S"));
		SAMEvidenceSource source = new SAMEvidenceSource(context, input, 0);
		ExtractEvidence e = new ExtractEvidence(context, source); e.process(EnumSet.allOf(ProcessStep.class)); e.close();
		assertEquals(1, new PerChr(context).getSC(source, "polyACGT").size());
	}
	@Test
	public void oea_should_be_located_in_rp_bam_for_chr() {
		ProcessingContext context = getCommandlineContext(true);
		createInput(OEA(1, 1, "100M", true));
		SAMEvidenceSource source = new SAMEvidenceSource(context, input, 0);
		ExtractEvidence e = new ExtractEvidence(context, source); e.process(EnumSet.allOf(ProcessStep.class)); e.close();
		assertEquals(1, new PerChr(context).getRP(source, "polyACGT").size());
	}
	@Test
	public void dp_should_be_located_in_rp_bam_for_chr() {
		ProcessingContext context = getCommandlineContext(true);
		createInput(DP(1, 1, "100M", true, 2, 2, "100M", true));
		SAMEvidenceSource source = new SAMEvidenceSource(context, input, 0);
		ExtractEvidence e = new ExtractEvidence(context, source); e.process(EnumSet.allOf(ProcessStep.class)); e.close();
		assertEquals(1, new PerChr(context).getRP(source, "polyACGT").size());
		assertEquals(1, new PerChr(context).getRP(source, "random").size());
	}
	@Test
	public void concordant_read_should_not_be_located_in_rp_or_scbam() {
		createInput(Read(1, 1, "100M"));
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, 0);
		ExtractEvidence e = new ExtractEvidence(getCommandlineContext(), source); e.process(EnumSet.allOf(ProcessStep.class)); e.close();
		assertEquals(0, getRP(source).size());
		assertEquals(0, getSC(source).size());
	}
	@Test
	public void oea_mate_should_be_located_in_mate_bam_for_chr() {
		ProcessingContext context = getCommandlineContext(true);
		createInput(OEA(1, 1, "100M", true));
		SAMEvidenceSource source = new SAMEvidenceSource(context, input, 0);
		ExtractEvidence e = new ExtractEvidence(context, source); e.process(EnumSet.allOf(ProcessStep.class)); e.close();
		List<SAMRecord> rs = new PerChr(context).getMate(source, "polyACGT");
		assertEquals(1, rs.size());
		SAMRecord mate = rs.get(0);
		assertTrue(mate.getReadUnmappedFlag());
		assertEquals(1, (int)mate.getMateReferenceIndex());
		assertEquals(1, mate.getMateAlignmentStart());
	}
	@Test
	public void dp_mate_should_be_located_in_mate_bam_for_chr() {
		ProcessingContext context = getCommandlineContext(true);
		createInput(DP(1, 1, "100M", true, 2, 2, "100M", true));
		SAMEvidenceSource source = new SAMEvidenceSource(context, input, 0);
		ExtractEvidence e = new ExtractEvidence(context, source); e.process(EnumSet.allOf(ProcessStep.class)); e.close();
		assertEquals(1, new PerChr(context).getMate(source, "polyACGT").size());
		assertEquals(1, new PerChr(context).getMate(source, "random").size());
		SAMRecord mate = new PerChr(context).getRP(source, "polyACGT").get(0);
		assertFalse(mate.getReadUnmappedFlag());
		assertEquals(2, mate.getMateAlignmentStart());
		assertEquals(2, mate.getMateAlignmentStart());
	}
	@Test
	public void mate_should_be_sorted_by_mate_coordinate() {
		createInput(DP(1, 1, "100M", true, 2, 5, "100M", true),
		   DP(1, 2, "100M", true, 2, 4, "100M", true),
		   DP(1, 3, "100M", true, 2, 6, "100M", true),
		   OEA(1, 4, "100M", false));
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, 0);
		ExtractEvidence e = new ExtractEvidence(getCommandlineContext(), source); e.process(EnumSet.allOf(ProcessStep.class)); e.close();
		List<SAMRecord> rs =  getMate(source);
		// polyACGT
		assertEquals(7, rs.size());
		assertEquals(1, rs.get(0).getMateAlignmentStart());
		assertEquals(2, rs.get(1).getMateAlignmentStart());
		assertEquals(3, rs.get(2).getMateAlignmentStart());
		assertEquals(4, rs.get(3).getMateAlignmentStart());
		// random 
		assertEquals(4, rs.get(4).getMateAlignmentStart());
		assertEquals(5, rs.get(5).getMateAlignmentStart());
		assertEquals(6, rs.get(6).getMateAlignmentStart());
	}
	@Test
	public void should_create_metrics() {
		createInput(RP(1, 1, 100, 10));
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, 0);
		ExtractEvidence e = new ExtractEvidence(getCommandlineContext(), source); e.process(EnumSet.allOf(ProcessStep.class)); e.close();
		assertTrue(getCommandlineContext().getFileSystemContext().getIdsvMetrics(input).exists());
		assertTrue(getCommandlineContext().getFileSystemContext().getInsertSizeMetrics(input).exists());
	}
	@Test
	public void should_process_metrics() {
		createInput(RP(0, 1, 2, 1), RP(0, 1, 7, 5));
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, 0);
		source.completeSteps(EnumSet.of(ProcessStep.CALCULATE_METRICS));
		assertEquals(5, source.getMaxReadLength());
		// 12345678901234567890
		// ----> <----
		assertEquals(11, source.getMaxConcordantFragmentSize());
		//assertEquals(PairOrientation.FR, metrics.getPairOrientation());
	}
	@Test
	public void should_set_NM_tag() {
		createInput(withSequence(S(POLY_A).substring(0, 100), Read(0, 1, "50M50S")));
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, 0);
		ExtractEvidence e = new ExtractEvidence(getCommandlineContext(), source); e.process(EnumSet.allOf(ProcessStep.class)); e.close();
		assertEquals(0, (int)getSC(source).get(0).getIntegerAttribute("NM"));
	}
	@Test
	public void should_write_long_sc_to_fastq() {
		createInput(ValidSC());
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, 0);
		ExtractEvidence e = new ExtractEvidence(getCommandlineContext(), source); e.process(EnumSet.allOf(ProcessStep.class)); e.close();
		
		List<FastqRecord> fastq = getFastqRecords(source);
		assertEquals(2, fastq.size());
		assertTrue(fastq.get(0).getReadHeader().contains("SC1"));
		assertTrue(fastq.get(1).getReadHeader().contains("SC1"));
		assertEquals("AAAAAAAAAAAAAAAAAAAAAAAAA", fastq.get(0).getReadString());
		// +33 encoding of mapping quality 5
		assertEquals("&&&&&&&&&&&&&&&&&&&&&&&&&", fastq.get(0).getBaseQualityString());
	}
	@Test
	public void should_not_write_indel_fastq() {
		createInput(Read(0, 1, "50M50D50M"));
		SAMEvidenceSource source = new SAMEvidenceSource(getCommandlineContext(), input, 0);
		ExtractEvidence e = new ExtractEvidence(getCommandlineContext(), source); e.process(EnumSet.allOf(ProcessStep.class)); e.close();
		
		List<FastqRecord> fastq = getFastqRecords(source);
		assertEquals(0, fastq.size());
	}
	@Test
	public void realign_min_mapq_should_filter_sc() {
		createInput(ValidSC());
		ProcessingContext pc = getCommandlineContext();
		pc.getConfig().minMapq = 16;
		SAMEvidenceSource source = new SAMEvidenceSource(pc, input, 0);
		ExtractEvidence e = new ExtractEvidence(pc, source); e.process(EnumSet.allOf(ProcessStep.class)); e.close();
		
		assertEquals(0, getFastqRecords(source).size());
	}
	@Test
	public void realign_max_mapq_should_filter_sc() {
		createInput(ValidSC());
		ProcessingContext pc = getCommandlineContext();
		pc.getConfig().maxMapq = 14;
		SAMEvidenceSource source = new SAMEvidenceSource(pc, input, 0);
		ExtractEvidence e = new ExtractEvidence(pc, source); e.process(EnumSet.allOf(ProcessStep.class)); e.close();
		
		assertEquals(0, getFastqRecords(source).size());
	}
	@Test
	public void realign_long_sc_length_should_filter_sc() {
		createInput(ValidSC());
		ProcessingContext pc = getCommandlineContext();
		pc.getRealignmentParameters().minLength = 26;
		SAMEvidenceSource source = new SAMEvidenceSource(pc, input, 0);
		ExtractEvidence e = new ExtractEvidence(pc, source); e.process(EnumSet.allOf(ProcessStep.class)); e.close();
		
		assertEquals(0, getFastqRecords(source).size());
	}
	@Test
	public void realign_min_id_should_filter_sc() {
		SAMRecord r = ValidSC();
		byte[] readBases = r.getReadBases();
		for (int i = 25; i < 51; i++) {
			// change just over half the bases from A to T
			readBases[i] = 'T';
		}
		r.setReadBases(readBases);
		createInput(RP(0, 1, 2, 1));
		ProcessingContext pc = getCommandlineContext();
		pc.getSoftClipParameters().minAnchorIdentity = 0.5f;
		SAMEvidenceSource source = new SAMEvidenceSource(pc, input, 0);
		ExtractEvidence e = new ExtractEvidence(pc, source); e.process(EnumSet.allOf(ProcessStep.class)); e.close();
		
		assertEquals(0, getFastqRecords(source).size());
	}
	@Test
	public void realign_min_qual_should_filter_sc() {
		createInput(ValidSC());
		ProcessingContext pc = getCommandlineContext();
		pc.getRealignmentParameters().minAverageQual = 6;
		SAMEvidenceSource source = new SAMEvidenceSource(pc, input, 0);
		ExtractEvidence e = new ExtractEvidence(pc, source); e.process(EnumSet.allOf(ProcessStep.class)); e.close();
		
		assertEquals(0, getFastqRecords(source).size());
	}
	@Test
	public void short_sc_should_not_be_considered_evidence() {
		createInput(Read(1, 1, "50M50S"));
		ProcessingContext pc = getCommandlineContext();
		pc.getSoftClipParameters().minLength = 51;
		SAMEvidenceSource source = new SAMEvidenceSource(pc, input, 0);
		ExtractEvidence e = new ExtractEvidence(pc, source); e.process(EnumSet.allOf(ProcessStep.class)); e.close();
		assertEquals(0, getSC(source).size());
	}
	@Test
	public void min_qual_should_only_consider_sc_qualities() {
		SAMRecord r = ValidSC();
		byte[] qual = r.getBaseQualities();
		for (int i = 25; i < 75; i++) {
			qual[i] = 40;
		}
		r.setBaseQualities(qual);
		createInput(r);
		ProcessingContext pc = getCommandlineContext();
		pc.getRealignmentParameters().minAverageQual = 6;
		SAMEvidenceSource source = new SAMEvidenceSource(pc, input, 0);
		ExtractEvidence e = new ExtractEvidence(pc, source); e.process(EnumSet.allOf(ProcessStep.class)); e.close();
		
		assertEquals(0, getFastqRecords(source).size());
	}
	@Test
	public void DREAM_unmapped_softclip_should_be_ignored() {
		SAMRecord r = Read(1, 1, "50M50S"); 
		r.setReadUnmappedFlag(true);
		createInput(r);
		ProcessingContext pc = getCommandlineContext();
		SAMEvidenceSource source = new SAMEvidenceSource(pc, input, 0);
		ExtractEvidence e = new ExtractEvidence(pc, source); e.process(ProcessStep.ALL_STEPS); e.close();
		assertEquals(0, getSC(source).size());
	}
	@Test
	public void should_extract_indels() {
		ProcessingContext context = getCommandlineContext(true);
		createInput(Read(1, 1, "50M50D50M"));
		SAMEvidenceSource source = new SAMEvidenceSource(context, input, 0);
		ExtractEvidence e = new ExtractEvidence(context, source); e.process(EnumSet.allOf(ProcessStep.class)); e.close();
		assertEquals(1, new PerChr(context).getSC(source, "polyACGT").size());
	}
	@Test
	public void should_not_extract_below_mapq_threshold() {
		createInput(withMapq(10, Read(1, 1, "50M50D50M")));
		ProcessingContext pc = getCommandlineContext(true);
		pc.getConfig().minMapq = 11;
		SAMEvidenceSource source = new SAMEvidenceSource(pc, input, 0);
		ExtractEvidence e = new ExtractEvidence(pc, source); e.process(EnumSet.allOf(ProcessStep.class)); e.close();
		assertEquals(0, new PerChr(pc).getSC(source, "polyACGT").size());
	}
}
