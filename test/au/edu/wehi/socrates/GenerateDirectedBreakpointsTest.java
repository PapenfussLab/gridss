package au.edu.wehi.socrates;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.List;

import net.sf.picard.cmdline.Option;
import net.sf.samtools.BAMRecordCodec;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordCoordinateComparator;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.util.SortingCollection;

import org.junit.After;
import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import com.google.common.collect.Lists;

public class GenerateDirectedBreakpointsTest extends CommandLineTest {
	public SAMRecord ValidSC() {
		SAMRecord r = Read(0, 1, "25S50M25S");
		r.setReadName("SC1");
		byte[] qual = new byte[100];
		for (int i = 0; i < qual.length; i++) {
			qual[i] = 5;
		}
		r.setBaseQualities(qual);
		r.setMappingQuality(5);
		return r;
	}
	@Test
	public void should_write_fastq() {
		createInput(RP(0, 1, 2, 1));
		extractEvidence();
		generateDirectedBreakpoints("polyA");
		assertTrue(fastq.exists());
	}
	@Test
	public void should_write_long_sc_to_fastq() {
		createInput(ValidSC());
		extractEvidence();
		generateDirectedBreakpoints("polyA");
		assertTrue(fastq.exists());
		assertEquals(2, getFastqRecords().size());
		assertEquals("fSC1", getFastqRecords().get(0).getReadHeader());
		assertEquals("bSC1", getFastqRecords().get(1).getReadHeader());
		assertEquals("AAAAAAAAAAAAAAAAAAAAAAAAA", getFastqRecords().get(0).getReadString());
		// +33 encoding of mapping quality 5
		assertEquals("&&&&&&&&&&&&&&&&&&&&&&&&&", getFastqRecords().get(0).getBaseQualityString());
	}
	@Test
	public void min_mapq_should_filter_sc() {
		createInput(ValidSC());
		extractEvidence();
		generateDirectedBreakpoints("polyA", 6, null, null, null, null);
		assertTrue(fastq.exists());
		assertEquals(0, getFastqRecords().size());
	}
	@Test
	public void long_sc_length_should_filter_sc() {
		createInput(ValidSC());
		extractEvidence();
		generateDirectedBreakpoints("polyA", null, 26, null, null, null);
		assertTrue(fastq.exists());
		assertEquals(0, getFastqRecords().size());
	}
	@Test
	public void min_id_should_filter_sc() {
		SAMRecord r = ValidSC();
		byte[] readBases = r.getReadBases();
		for (int i = 25; i < 51; i++) {
			// change just over half the bases from A to T
			readBases[i] = 'T';
		}
		r.setReadBases(readBases);
		createInput(RP(0, 1, 2, 1));
		extractEvidence();
		generateDirectedBreakpoints("polyA", null, null, 50f, null, null);
		assertTrue(fastq.exists());
		assertEquals(0, getFastqRecords().size());
	}
	@Test
	public void min_qual_should_filter_sc() {
		createInput(ValidSC());
		extractEvidence();
		generateDirectedBreakpoints("polyA", null, null, null, 6f, null);
		assertTrue(fastq.exists());
		assertEquals(0, getFastqRecords().size());
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
		extractEvidence();
		generateDirectedBreakpoints("polyA", null, null, null, 6f, null);
		assertTrue(fastq.exists());
		assertEquals(0, getFastqRecords().size());
	}
	@Test
	public void debruijn_should_generate_fastq() {
		createInput();
		extractEvidence();
		generateDirectedBreakpoints("polyA", null, null, null, null, null);
		assertTrue(fastq.exists());
		assertEquals(1, getFastqRecords().size());
		assertTrue("TODO: de bruijn graph", false);
	}
	@Test
	public void debruijn_should_generate_vcf() {
		createInput(RP(0, 1, 2, 1));
		extractEvidence();
		generateDirectedBreakpoints("polyA", null, null, null, null, null);
		assertEquals(1, getVcf(".socrates.polyA.breakend.vcf").size());
	}
	@Test
	public void should_write_vcf() {
		createInput(RP(0, 1, 2, 1));
		extractEvidence();
		generateDirectedBreakpoints("polyA");
		shouldExist("input.bam.socrates.polyA.breakend.vcf");
	}
}
