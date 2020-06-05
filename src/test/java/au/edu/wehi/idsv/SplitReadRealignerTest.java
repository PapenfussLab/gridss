package au.edu.wehi.idsv;

import au.edu.wehi.idsv.alignment.ExternalAlignerTests;
import au.edu.wehi.idsv.alignment.ExternalProcessStreamingAligner;
import au.edu.wehi.idsv.alignment.FastqAligner;
import au.edu.wehi.idsv.picard.BufferedReferenceSequenceFile;
import au.edu.wehi.idsv.picard.SynchronousReferenceLookupAdapter;
import au.edu.wehi.idsv.util.UngroupingIterator;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import joptsimple.internal.Strings;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.junit.experimental.categories.Category;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;

import static org.junit.Assert.*;

public abstract class SplitReadRealignerTest extends IntermediateFilesTest {
	private SplitReadRealigner getSplitReadFastqExtractor(
			int minSoftClipLength,
			float minClipQuality,
			boolean processSecondaryAlignments,
			boolean realignExistingSplitReads,
			boolean realignEntireRecord,
			boolean realignAnchoringBases,
			EvidenceIdentifierGenerator eidgen) {
		SplitReadRealigner srr = new StreamingSplitReadRealigner(getContext(), null, 1);
		srr.setMinSoftClipLength(minSoftClipLength);
		srr.setMinSoftClipQuality(minClipQuality);
		srr.setProcessSecondaryAlignments(processSecondaryAlignments);
		srr.setRealignExistingSplitReads(realignExistingSplitReads);
		srr.setRealignEntireRecord(realignEntireRecord);
		srr.setRealignAnchoringBases(realignAnchoringBases);
		srr.setEvidenceIdentifierGenerator(eidgen);
		return srr;
	}
	private Iterator<FastqRecord> getSplitReadFastqExtractionIterator(
			Iterator<SAMRecord> it,
			boolean recursive,
			int minSoftClipLength,
			float minClipQuality,
			boolean processSecondaryAlignments,
			boolean realignExistingSplitReads,
			boolean realignEntireRecord,
			boolean realignAnchoringBases,
			EvidenceIdentifierGenerator eidgen) {
		SplitReadRealigner srfe = getSplitReadFastqExtractor(minSoftClipLength, minClipQuality, processSecondaryAlignments, realignExistingSplitReads, realignEntireRecord, realignAnchoringBases, eidgen);
		return new UngroupingIterator<>(Iterators.transform(it, (SAMRecord r) -> srfe.extract(r, recursive)));
	}
	@Test
	public void should_extract_full_sequence() {
		SplitReadRealigner srfe = getSplitReadFastqExtractor(1, 0, false, false, true, false, new StringEvidenceIdentifierGenerator());
		List<FastqRecord> result = srfe.extract(Read(0, 1, "25M75S"), false);
		Assert.assertEquals(1, result.size());
		Assert.assertEquals(100, result.get(0).getReadLength());
	}
	@Test
	public void should_extract_remaining_sequence() {
		SplitReadRealigner srfe  = getSplitReadFastqExtractor(1, 0, false, false, true, false, new StringEvidenceIdentifierGenerator());
		FastqRecord fqr = srfe.extract(Read(0, 1, "25M75S"), false).get(0);
		SAMRecord aligned = new SAMRecord(getHeader());
		aligned.setReadBases(fqr.getReadBases());
		aligned.setReadName(fqr.getReadName());
		aligned.setBaseQualities(fqr.getBaseQualities());
		aligned.setReferenceIndex(0);
		aligned.setAlignmentStart(0);
		aligned.setCigarString("10S40M50S");
		srfe = getSplitReadFastqExtractor(1, 0, false, false, true, false, new StringEvidenceIdentifierGenerator());
		List<FastqRecord> result = srfe.extract(aligned, true);
		Assert.assertEquals(2, result.size());
		Assert.assertEquals(10, result.get(0).getReadLength());
		Assert.assertEquals(50, result.get(1).getReadLength());
	}
	@Test
	public void should_pad_anchor_with_N() {
		SplitReadRealigner srfe  = getSplitReadFastqExtractor( 1, 0, false, false, false, true, new StringEvidenceIdentifierGenerator());
		List<FastqRecord> result = srfe.extract(Read(0, 1, "1S2M3S"), false);
		Assert.assertEquals(3, result.size());
		Assert.assertEquals("NAANNN", result.get(2).getReadString());
	}
	@Test
	public void should_extract_above_min_soft_clip_length() {
		assertEquals(0, Iterators.size(getSplitReadFastqExtractionIterator(ImmutableList.of(Read(0, 1, "1S10M")).iterator(), false, 2, 0, true, false, false, false, getContext().getEvidenceIDGenerator())));
		assertEquals(1, Iterators.size(getSplitReadFastqExtractionIterator(ImmutableList.of(Read(0, 1, "2S10M")).iterator(), false, 2, 0, true, false, false, false, getContext().getEvidenceIDGenerator())));
	}

	@Test
	public void should_filter_secondary() {
		SAMRecord r = Read(0, 1, "10S10M");
		r.setSecondaryAlignment(true);
		assertEquals(0, Iterators.size(getSplitReadFastqExtractionIterator(ImmutableList.of(r).iterator(), false, 1, 0, false, false, false, false, getContext().getEvidenceIDGenerator())));
		assertEquals(1, Iterators.size(getSplitReadFastqExtractionIterator(ImmutableList.of(r).iterator(), false, 1, 0, true, false, false, false, getContext().getEvidenceIDGenerator())));
	}
	@Test
	public void should_filter_supplemenary() {
		SAMRecord r = Read(0, 1, "10S10M");
		r.setSupplementaryAlignmentFlag(true);
		assertEquals(0, Iterators.size(getSplitReadFastqExtractionIterator(ImmutableList.of(r).iterator(), false, 1, 0, true, false, false, false, getContext().getEvidenceIDGenerator())));
	}
	@Test
	public void should_filter_unmapped() {
		SAMRecord r = Read(0, 1, "10S10M");
		r.setReadUnmappedFlag(true);
		assertEquals(0, Iterators.size(getSplitReadFastqExtractionIterator(ImmutableList.of(r).iterator(), false, 1, 0, true, false, false, false, getContext().getEvidenceIDGenerator())));
	}
	@Test
	public void should_filter_existing_split_reads_with_SA_tag() {
		SAMRecord r = Read(0, 1, "10S10M");
		r.setAttribute("SA", "polyA,100,+,10M10S,0,0");
		assertEquals(0, Iterators.size(getSplitReadFastqExtractionIterator(ImmutableList.of(r).iterator(), false, 1, 0, true, false, false, false, getContext().getEvidenceIDGenerator())));
	}
	@Test
	public void should_not_filter_existing_split_reads_with_SA_tag_if_realignExistingSplitReads_is_true() {
		SAMRecord r = Read(0, 1, "10S10M");
		r.setAttribute("SA", "polyA,100,+,10M10S,0,0");
		assertEquals(1, Iterators.size(getSplitReadFastqExtractionIterator(ImmutableList.of(r).iterator(), false, 1, 0, true, true, false, false, getContext().getEvidenceIDGenerator())));
	}
	@Test
	public void should_filter_under_average_mapq() {
		SAMRecord r = Read(0, 1, "2S1M");
		r.setBaseQualities(new byte[] { 0, 10, 40 } );
		assertEquals(1, Iterators.size(getSplitReadFastqExtractionIterator(ImmutableList.of(r).iterator(), false, 1, 4, true, false, false, false, getContext().getEvidenceIDGenerator())));
		assertEquals(1, Iterators.size(getSplitReadFastqExtractionIterator(ImmutableList.of(r).iterator(), false, 1, 5, true, false, false, false, getContext().getEvidenceIDGenerator())));
		assertEquals(0, Iterators.size(getSplitReadFastqExtractionIterator(ImmutableList.of(r).iterator(), false, 1, 6, true, false, false, false, getContext().getEvidenceIDGenerator())));
	}

	protected SplitReadRealigner srr = null;
	@Before
	@Override
	public void setup() throws IOException {
		super.setup();
		output = testFolder.newFile("out.sam");
		output.delete();
		srr = createAligner();
	}
	protected abstract SplitReadRealigner createAligner();
	@Test
	public void should_output_new_records_to_modified_file() throws IOException {
		File outputModified = new File(output.toPath() + ".modified.bam");
		createInput(withSequence(
				"CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACC" +
					 "GAACAGGCTTTTTGTGCCATGCCTGATGCCATACTCCGAACACATGCCTA",
				Read(2, 1, "50M50S")));
		srr.setRealignEntireRecord(true);
		srr.createSupplementaryAlignments(input, output, outputModified);
		assertEquals(1, getRecords(output).size());
		assertEquals(1, getRecords(outputModified).size());
	}
	@Test
	public void should_output_moved_records_to_different_file() throws IOException {
		File outputModified = new File(output.toPath() + ".modified.bam");
		createInput(withSequence(
				"CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACC" +
					 "GAACAGGCTTTTTGTGCCATGCCTGATGCCATACTCCGAACACATGCCTA",
				Read(2, 2, "50M50S")));
		srr.setRealignEntireRecord(true);
		srr.createSupplementaryAlignments(input, output, outputModified);
		List<SAMRecord> result = getRecords(outputModified);
		assertEquals(0, getRecords(output).size());
		assertEquals(2, getRecords(outputModified).size());
	}
	@Test
	public void should_not_call_aligner_unless_required() throws IOException {
		createInput(Read(0, 1, "50M"));
		FastqAligner aligner = null;
		srr.createSupplementaryAlignments(input, output, output);
		assertEquals(1, getRecords(output).size());
	}
	@Test
	public void should_output_sorted_merge() throws IOException {
		SAMRecord r = Read(2, 1, "50S50M");
		r.setReadBases(B(S(RANDOM).substring(100, 150) + S(RANDOM).substring(0, 50)));
		r.setReadName("r");
		
		createBAM(input, SortOrder.coordinate, r);
		srr.createSupplementaryAlignments(input, output, output);
		
		List<SAMRecord> result = getRecords(output);
		assertEquals(2, result.size());
		assertEquals(1, result.get(0).getAlignmentStart());
		assertEquals(101, result.get(1).getAlignmentStart());
	}
	@Test
	public void should_recusively_align() throws IOException {
		SAMRecord r = Read(2, 1, "50S50M");
		r.setReadBases(B(S(RANDOM).substring(125, 150) + S(RANDOM).substring(75, 100) + S(RANDOM).substring(0, 50)));
		r.setReadName("r");
		
		createBAM(input, SortOrder.coordinate, r);
		srr.createSupplementaryAlignments(input, output, output);
		
		List<SAMRecord> result = getRecords(output);
		assertEquals(3, result.size());
		assertEquals(1, result.get(0).getAlignmentStart());
		assertEquals(76, result.get(1).getAlignmentStart());
		assertEquals(126, result.get(2).getAlignmentStart());
	}
	@Test
	public void should_match_input_sort_order() throws IOException {
		SAMRecord r = Read(2, 1, "50S50M");
		r.setReadBases(B(S(RANDOM).substring(125, 150) + S(RANDOM).substring(75, 100) + S(RANDOM).substring(0, 50)));
		r.setReadName("r");
		List<SAMRecord> result;

		createBAM(input, SortOrder.coordinate, r);
		srr.createSupplementaryAlignments(input, output, output);
		result = getRecords(output);
		assertEquals(3, result.size());
		assertTrue(Ordering.from(SortOrder.coordinate.getComparatorInstance()).isOrdered(result));
		
		createBAM(input, SortOrder.queryname, r);
		srr.createSupplementaryAlignments(input, output, output);
		result = getRecords(output);
		assertEquals(3, result.size());
		assertTrue(Ordering.from(SortOrder.queryname.getComparatorInstance()).isOrdered(result));
		
		createBAM(input, SortOrder.unsorted, r);
		srr.createSupplementaryAlignments(input, output, output);
		result = getRecords(output);
		assertEquals(3, result.size());
	}
	@Test
	public void fastq_keys_should_be_unique() throws IOException {
		SAMRecord r = Read(0, 1, "1S1M1S");
		List<FastqRecord> fq = srr.extract(r, false);
		assertEquals(2, fq.size());
		assertNotEquals(fq.get(0).getReadName(), fq.get(1).getReadName());
	}
	@Test
	public void fastq_keys_should_not_exceed_254_character_BAM_limit() throws IOException {
		SAMRecord r = Read(0, 1, "1S1M1S");
		r.setReadName(S("R", 254));
		List<FastqRecord> fq = srr.extract(r, false);
		assertEquals(2, fq.size());
		for (FastqRecord fqr : fq) {
			assertTrue(fqr.getReadName().length() <= 254);
		}
	}
	@Test
	public void should_realign_multiple_times() throws IOException {
		SAMRecord r = Read(2, 1, "50M150S");
		r.setReadBases(B(S(RANDOM).substring(0, 50) + S(RANDOM).substring(100, 150) + S(RANDOM).substring(200, 250) + S(RANDOM).substring(300, 350)));
		createBAM(input, getHeader(), r);

		srr.createSupplementaryAlignments(input, output, output);
		List<SAMRecord> list = getRecords(output);
		assertEquals(4, list.size());
	}

	@Test
	public void realign_entire_read_should_drop_existing_supplementary_alignments() throws IOException {
		SAMRecord r1 = Read(2, 1, "50M100S");
		SAMRecord r1s = Read(0, 1, "50S100M");
		r1.setReadBases(B(S(RANDOM).substring(0, 50) + S(RANDOM).substring(100, 150) + S(RANDOM).substring(200, 250)));
		r1s.setReadBases(r1.getReadBases());
		r1s.setSupplementaryAlignmentFlag(true);
		r1s.setReadName(r1.getReadName());
		createBAM(input, getHeader(), r1, r1s);

		srr.setRealignEntireRecord(true);
		assertFalse(srr.shouldDropInputRecord(r1));
		assertTrue(srr.shouldDropInputRecord(r1s));
	}

	@Test
	public void realign_existing_split_reads_should_drop_existing_supplementary_alignments() throws IOException {
		SAMRecord r1 = Read(2, 1, "50M100S");
		SAMRecord r1s = Read(0, 1, "50S100M");
		r1.setReadBases(B(S(RANDOM).substring(0, 50) + S(RANDOM).substring(100, 150) + S(RANDOM).substring(200, 250)));
		r1s.setReadBases(r1.getReadBases());
		r1s.setSupplementaryAlignmentFlag(true);
		r1s.setReadName(r1.getReadName());
		createBAM(input, getHeader(), r1, r1s);

		srr.setRealignExistingSplitReads(true);
		assertFalse(srr.shouldDropInputRecord(r1));
		assertTrue(srr.shouldDropInputRecord(r1s));
	}


	@Test
	public void realign_full_read_should_drop_existing_supp_records() throws IOException {
		SAMRecord r = new SAMRecord(getHeader());
		r.setReferenceIndex(0);
		r.setAlignmentStart(1);
		r.setCigarString("30S30M");
		r.setReadBases(B("AAAAAAAAAATTTTTTTTTTGGGGGGGGGGAAAAAAAAAATTTTTTTTTTGGGGGGGGGG"));
		r.setBaseQualities(B("AAAAAAAAAATTTTTTTTTTGGGGGGGGGGAAAAAAAAAATTTTTTTTTTGGGGGGGGGG"));
		r.setReadName("supp");
		r.setSupplementaryAlignmentFlag(true);

		createBAM(input, getHeader(), r);
		SplitReadRealigner aligner = createAligner();
		srr.setRealignEntireRecord(true);
		srr.createSupplementaryAlignments(input, output, output);
		List<SAMRecord> list = getRecords(output);
		assertEquals(0, list.size());
	}

	@Test
	public void realign_anchoring_bases_should_only_write_primary_or_overlapping_anchoring_base_alignment() throws IOException {
		SAMRecord r = new SAMRecord(getHeader());
		r.setReferenceIndex(2);
		r.setAlignmentStart(100);
		r.setCigarString("100S300M");
		r.setReadBases(B(S(RANDOM).substring(300, 400) + S(RANDOM).substring(100, 150) + S(RANDOM).substring(1000, 1250)));
		r.setBaseQualities(getPolyA(200));
		r.setReadName("anchor_realign_test");

		createBAM(input, getHeader(), r);
		SplitReadRealigner aligner = createAligner();
		srr.setRealignEntireRecord(true);
		srr.createSupplementaryAlignments(input, output, output);
		List<SAMRecord> list = getRecords(output);
		assertEquals(2, list.size());
		assertEquals(false, list.get(0).getSupplementaryAlignmentFlag());
		assertEquals("100S50M250S", list.get(0).getCigarString());
		assertEquals("100S300M", list.get(0).getStringAttribute("OA"));
		assertEquals(true, list.get(1).getSupplementaryAlignmentFlag());
		assertEquals("100M300S", list.get(1).getCigarString());
		assertEquals("100S300M", list.get(1).getStringAttribute("OA"));
	}

	/**
	 * bwa mem primary alignment record over-aligns across apparent SNVs and an indels which should just be placed on the other side  
	 * @throws IOException
	 */
	@Test
	@Category({Hg19Tests.class, ExternalAlignerTests.class})
	public void realign_should_not_overalign_one_side() throws IOException {
		File hg19 = Hg19Tests.findHg19Reference();
		SynchronousReferenceLookupAdapter ref = new SynchronousReferenceLookupAdapter(new IndexedFastaSequenceFile(hg19));
		ExternalProcessStreamingAligner aligner = new ExternalProcessStreamingAligner(SamReaderFactory.makeDefault(), ExternalAlignerTests.COMMAND_LINE, hg19, 4, new IndexedFastaSequenceFile(ExternalAlignerTests.REFERENCE).getSequenceDictionary());
		ProcessingContext pc = new ProcessingContext(new FileSystemContext(testFolder.getRoot(), 500000), hg19, ref, Lists.newArrayList(), getConfig(testFolder.getRoot()));
		SplitReadRealigner srr = new StreamingSplitReadRealigner(pc, aligner, 1);
		srr.setRealignEntireRecord(true);
		SAMFileHeader header = new SAMFileHeader();
		header.setSequenceDictionary(ref.getSequenceDictionary());
		header.setSortOrder(SortOrder.coordinate);
		
		SAMRecord r = new SAMRecord(header);
		r.setReferenceIndex(0);
		r.setAlignmentStart(1);
		r.setCigarString("800S100M");
		r.setReadBases(    B("TTGTTATTTGTAATAGCATCTAACTGAAATCTAGCATTTAGAAATCTCTCTATCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCCTTCTTTCTTTCTTTCTTCTCTCTCTCTCTCTCTTTCTTTTGATGGAATTTTGCTCTTGTCGCCCAGACTGGAGTGCAATGGTGCGATCTTGGCTCACAGCAACCTCTGCCTCCCAGGTTCAAGCAATTCTTCTGTCTCAGCCTCCCAAGGAGCTGGGATTACAGGCACATGCCACCATGCCAGGCTAATTTTTGTATTTTTAGTAGAGACGGGGTTTCACCATGTTGGCCAAGCTGGTCTTGAACTCAAGCTCATGTTCCACCTGCCTTGGCCTCTCTACGTGCTGGGATTACAGGCATGAGCTACAGCGCTCAGCCGAAAGATTAATTTTTTAAAAATGCCTCCTGTGACCTCCATTTCTCTTCTACTTAGAGTCAGAGAGTGTTAAAACTTACAGACCAACTTAATCCAAGCCCTCCATTGGGCAGATTGGGATGCTGAGGAGAAGATGGAGGGGGATTTGTTTTTGGTAGCCTAGGGGTCCCACAGGAGAAAAGAAAGGAGGCAGTAGGACAGGAAAGGAAAAATTTGGAGGCAGGTTTGGAGGGTAAGGGGTGAAATAGCTGGCTTGGGCTACCAGAGGATCTGGAGCAGGCACCCTTGGCCTTTATAGTTGCAGGGACCCATTTGGCAGGCTGGTGAAGTCTATGAATCCCTTCTCAGGAAGTGTTTTAAATGCATGAAATAAAATGAAAGACATTGGATTCTAAAGAAAACCAAATATATTGAAATATAGTTACCAAAATATTTTAAAAACAAATGTGTGATTTAGTAATACAAGTACTTCTTTATTAATGCATT"));
		r.setBaseQualities(B("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
		r.setReadName("read");
		
		createBAM(input, header, r);
		srr.createSupplementaryAlignments(input, output, output);
		List<SAMRecord> list = getRecords(output);
		assertEquals(2, list.size());
		Assert.assertEquals("282S618M", list.get(0).getCigarString());
	}

}
