package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;

import au.edu.wehi.idsv.configuration.GridssConfiguration;
import com.google.common.collect.ImmutableList;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;
import org.junit.experimental.categories.Category;

import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;

import au.edu.wehi.idsv.alignment.AlignerFactory;
import au.edu.wehi.idsv.alignment.ExternalAlignerTests;
import au.edu.wehi.idsv.alignment.ExternalProcessStreamingAligner;
import au.edu.wehi.idsv.alignment.FastqAligner;
import au.edu.wehi.idsv.alignment.SmithWatermanFastqAligner;
import au.edu.wehi.idsv.alignment.StreamingAligner;
import au.edu.wehi.idsv.picard.BufferedReferenceSequenceFile;
import au.edu.wehi.idsv.picard.SynchronousReferenceLookupAdapter;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;

public class SplitReadRealignerTest extends IntermediateFilesTest {
	private static final SmithWatermanFastqAligner aligner = new SmithWatermanFastqAligner(AlignerFactory.create(), 2);
	@Before
	@Override
	public void setup() throws IOException {
		super.setup();
		output = testFolder.newFile("out.sam");
		output.delete();
	}
	@Test
	public void should_output_new_records_to_modified_file() throws IOException {
		File outputModified = new File(output.toPath() + ".modified.bam");
		createInput(withSequence(
				"CATTAATCGCAAGAGCGGGTTGTATTCGACGCCAAGTCAGCTGAAGCACC" +
					 "GAACAGGCTTTTTGTGCCATGCCTGATGCCATACTCCGAACACATGCCTA",
				Read(2, 1, "50M50S")));
		SplitReadRealigner srr = new SplitReadRealigner(getContext());
		srr.setRealignEntireRecord(true);
		srr.createSupplementaryAlignments(aligner, input, output, outputModified, true);
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
		SplitReadRealigner srr = new SplitReadRealigner(getContext());
		srr.setRealignEntireRecord(true);
		srr.createSupplementaryAlignments(aligner, input, output, outputModified, true);
		List<SAMRecord> result = getRecords(outputModified);
		assertEquals(0, getRecords(output).size());
		assertEquals(2, getRecords(outputModified).size());
	}
	@Test
	public void should_not_call_aligner_unless_required() throws IOException {
		createInput(Read(0, 1, "50M"));
		SplitReadRealigner srr = new SplitReadRealigner(getContext());
		FastqAligner aligner = null;
		srr.createSupplementaryAlignments(aligner, input, output, output, true);
		assertEquals(1, getRecords(output).size());
	}
	@Test
	public void should_output_sorted_merge() throws IOException {
		SAMRecord r = Read(2, 1, "50S50M");
		r.setReadBases(B(S(RANDOM).substring(100, 150) + S(RANDOM).substring(0, 50)));
		r.setReadName("r");
		
		createBAM(input, SortOrder.coordinate, r);
		SplitReadRealigner srr = new SplitReadRealigner(getContext());
		srr.createSupplementaryAlignments(aligner, input, output, output, true);
		
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
		SplitReadRealigner srr = new SplitReadRealigner(getContext());
		srr.createSupplementaryAlignments(aligner, input, output, output, true);
		
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
		SplitReadRealigner srr;
		
		createBAM(input, SortOrder.coordinate, r);
		srr = new SplitReadRealigner(getContext());
		srr.createSupplementaryAlignments(aligner, input, output, output, true);
		result = getRecords(output);
		assertEquals(3, result.size());
		assertTrue(Ordering.from(SortOrder.coordinate.getComparatorInstance()).isOrdered(result));
		
		createBAM(input, SortOrder.queryname, r);
		srr = new SplitReadRealigner(getContext());
		srr.createSupplementaryAlignments(aligner, input, output, output, true);
		result = getRecords(output);
		assertEquals(3, result.size());
		assertTrue(Ordering.from(SortOrder.queryname.getComparatorInstance()).isOrdered(result));
		
		createBAM(input, SortOrder.unsorted, r);
		srr = new SplitReadRealigner(getContext());
		srr.createSupplementaryAlignments(aligner, input, output, output, true);
		result = getRecords(output);
		assertEquals(3, result.size());
	}
	@Test
	public void fastq_keys_should_be_unique() throws IOException {
		createBAM(input, SortOrder.coordinate, Read(0, 1, "1S1M1S"));
		File fq1 = getContext().getFileSystemContext().getRealignmentFastq(input, 0);
		new SplitReadRealigner(getContext()).createSupplementaryAlignmentFastq(input, fq1, false);
		List<FastqRecord> fq = getFastqRecords(fq1);
		assertEquals(2, fq.size());
		assertNotEquals(fq.get(0).getReadName(), fq.get(1).getReadName());
	}
	@Test
	public void fastq_keys_should_not_exceed_254_character_BAM_limit() throws IOException {
		SAMRecord r = Read(0, 1, "1S1M1S");
		r.setReadName(S("R", 254));
		createBAM(input, SortOrder.coordinate, r);
		File fq1 = getContext().getFileSystemContext().getRealignmentFastq(input, 0);
		SplitReadRealigner srr = new SplitReadRealigner(getContext());
		srr.createSupplementaryAlignmentFastq(input, fq1, false);
		List<FastqRecord> fq = getFastqRecords(fq1);
		for (FastqRecord fqr : fq) {
			assertTrue(fqr.getReadName().length() <= 254);
		}
	}
	public class StubStreamingAligner implements StreamingAligner {
		private int out = 0;
		private int in = 0;
		private SAMRecord[] alignments;
		public StubStreamingAligner(SAMRecord... alignments) {
			this.alignments = alignments;
		}
		@Override
		public void asyncAlign(FastqRecord fq) throws IOException {
			this.alignments[in++].setReadName(fq.getReadName());
		}
		@Override
		public void flush() throws IOException { }

		@Override
		public boolean hasAlignmentRecord() {
			return in > out & out < alignments.length;
		}

		@Override
		public SAMRecord getAlignment() {
			return alignments[out++];
		}

		@Override
		public void close() throws IOException { }
		@Override
		public int processedAlignmentRecords() {
			// alignment is immediate
			if (!hasAlignmentRecord()) return 0;
			return out - in;
		}
		@Override
		public int outstandingAlignmentRecord() {
			return 0;
		}
	}
	@Test
	public void streaming_should_realign_multiple_times() throws IOException {
		SAMRecord r0 = Read(0, 100, "10M30S");
		r0.setReadName("r");
		SAMRecord r1 = Read(0, 200, "10M20S");
		r1.setAttribute("SA", "this record is a split read realignment");
		SAMRecord r2 = Read(0, 200, "10M10S");
		
		SAMRecord r3 = Read(0, 200, "10M");
		
		createBAM(input, SortOrder.coordinate, r0);
		SplitReadRealigner srr = new SplitReadRealigner(getContext());
		srr.createSupplementaryAlignments(new StubStreamingAligner(r1, r2, r3), input, output, output, true, 100);
		List<SAMRecord> list = getRecords(output);
		assertEquals(4, list.size());
	}
	@Test
	@Category(ExternalAlignerTests.class)
	@Ignore("Working 2018-04-08. Currently need to manual check # restarts of external aligner to actually test this functionality. Needs a delayed return stub to test properly.")
	public void streaming_should_limit_outstanding_records_to_buffer_size() throws IOException, CloneNotSupportedException {
		ExternalProcessStreamingAligner aligner = new ExternalProcessStreamingAligner(SamReaderFactory.makeDefault(), ExternalAlignerTests.COMMAND_LINE, ExternalAlignerTests.REFERENCE, 4, new IndexedFastaSequenceFile(ExternalAlignerTests.REFERENCE).getSequenceDictionary());
		BufferedReferenceSequenceFile lookup = new BufferedReferenceSequenceFile(ReferenceSequenceFileFactory.getReferenceSequenceFile(ExternalAlignerTests.REFERENCE));
		ProcessingContext pc = new ProcessingContext(new FileSystemContext(testFolder.getRoot(), 500000), ExternalAlignerTests.REFERENCE, lookup, Lists.newArrayList(), getConfig(testFolder.getRoot()));
		SplitReadRealigner srr = new SplitReadRealigner(pc);
		
		SAMFileHeader header = new SAMFileHeader();
		header.setSequenceDictionary(lookup.getSequenceDictionary());
		header.setSortOrder(SortOrder.coordinate);
		
		
		SAMRecord r = new SAMRecord(header);
		r.setReferenceIndex(0);
		r.setAlignmentStart(1399998);
		r.setCigarString("301S103M");
		r.setReadBases(B("GGATATATAGGGATAGAAGCTTGAATAGTCTGGACATATATTTGTATTGAAATACAAATGTAAGATTTCAGTTAATCAATTTAAACATTTTTATTTTCAAGGGCTTCCAGCGTCCACTTCCTACGGCAAGCAGGAGGAGACAAGCGCCACCCTGCGCTCGCGGAGCCGACCCCGGCTCTCCCCTCCCGTGGCCGCAGGGGTCTGACAGAAAGGGGTCACTAATCTACTTGGCCTTTTGAGGACTGATCCTTAAGAATAATTTTTTTTTTTTTATGATCTTGAAGGCTGAGAAGTATTAGAGTAGGTTTTTTTCTCCTTCATAAGGCCAGATTCTTCTTTCTGTCACAGATTTCAAGTCCCCGCCTCAGCAGCCTTTCACTGTCAGTTCTTTCTCACGTGACCCT"));
		r.setBaseQualities(B("?????BBBB@DEDDDDGGGGGEIEHIHEFHIHIIEHHIEIIIIIIEHII?HHFHHHHDIHIHEHHFIIBCHI=GHIH@HFCEIGIHIDHHHGCIIHDHHFA?????BBBB@DEDDDDGGGGGEIEHIHEFHIHIIEHHIEIIIIIIEHII?HHFHHHHDIHIHEHHFIIBCHI=GHIH@HFCEIGIHIDHHHGCIIHDHHFA?????BBBB@DEDDDDGGGGGEIEHIHEFHIHIIEHHIEIIIIIIEHII?HHFHHHHDIHIHEHHFIIBCHI=GHIH@HFCEIGIHIDHHHGCIIHDHHFA?????BBBB@DEDDDDGGGGGEIEHIHEFHIHIIEHHIEIIIIIIEHII?HHFHHHHDIHIHEHHFIIBCHI=GHIH@HFCEIGIHIDHHHGCIIHDHHFA"));
		r.setReadName("four_way_split_read");
		
		SAMRecord[] arr = new SAMRecord[64];
		for (int i = 0; i < arr.length; i++) {
			arr[i] = (SAMRecord)r.clone();
			arr[i].setReadName(Integer.toString(i));
		}
		createBAM(input, header, arr);
		
		srr.createSupplementaryAlignments(aligner, input, output, output, true, 16);
		List<SAMRecord> list = getRecords(output);
		assertEquals(4, list.size());
	}
	@Test
	@Category(ExternalAlignerTests.class)
	public void should_realign_multiple_times() throws IOException {
		ExternalProcessStreamingAligner aligner = new ExternalProcessStreamingAligner(SamReaderFactory.makeDefault(), ExternalAlignerTests.COMMAND_LINE, ExternalAlignerTests.REFERENCE, 4, new IndexedFastaSequenceFile(ExternalAlignerTests.REFERENCE).getSequenceDictionary());
		BufferedReferenceSequenceFile lookup = new BufferedReferenceSequenceFile(ReferenceSequenceFileFactory.getReferenceSequenceFile(ExternalAlignerTests.REFERENCE));
		ProcessingContext pc = new ProcessingContext(new FileSystemContext(testFolder.getRoot(), 500000), ExternalAlignerTests.REFERENCE, lookup, Lists.newArrayList(), getConfig(testFolder.getRoot()));
		SplitReadRealigner srr = new SplitReadRealigner(pc);
		
		SAMFileHeader header = new SAMFileHeader();
		header.setSequenceDictionary(lookup.getSequenceDictionary());
		header.setSortOrder(SortOrder.coordinate);
		
		
		SAMRecord r = new SAMRecord(header);
		r.setReferenceIndex(0);
		r.setAlignmentStart(1399998);
		r.setCigarString("301S103M");
		r.setReadBases(B("GGATATATAGGGATAGAAGCTTGAATAGTCTGGACATATATTTGTATTGAAATACAAATGTAAGATTTCAGTTAATCAATTTAAACATTTTTATTTTCAAGGGCTTCCAGCGTCCACTTCCTACGGCAAGCAGGAGGAGACAAGCGCCACCCTGCGCTCGCGGAGCCGACCCCGGCTCTCCCCTCCCGTGGCCGCAGGGGTCTGACAGAAAGGGGTCACTAATCTACTTGGCCTTTTGAGGACTGATCCTTAAGAATAATTTTTTTTTTTTTATGATCTTGAAGGCTGAGAAGTATTAGAGTAGGTTTTTTTCTCCTTCATAAGGCCAGATTCTTCTTTCTGTCACAGATTTCAAGTCCCCGCCTCAGCAGCCTTTCACTGTCAGTTCTTTCTCACGTGACCCT"));
		r.setBaseQualities(B("?????BBBB@DEDDDDGGGGGEIEHIHEFHIHIIEHHIEIIIIIIEHII?HHFHHHHDIHIHEHHFIIBCHI=GHIH@HFCEIGIHIDHHHGCIIHDHHFA?????BBBB@DEDDDDGGGGGEIEHIHEFHIHIIEHHIEIIIIIIEHII?HHFHHHHDIHIHEHHFIIBCHI=GHIH@HFCEIGIHIDHHHGCIIHDHHFA?????BBBB@DEDDDDGGGGGEIEHIHEFHIHIIEHHIEIIIIIIEHII?HHFHHHHDIHIHEHHFIIBCHI=GHIH@HFCEIGIHIDHHHGCIIHDHHFA?????BBBB@DEDDDDGGGGGEIEHIHEFHIHIIEHHIEIIIIIIEHII?HHFHHHHDIHIHEHHFIIBCHI=GHIH@HFCEIGIHIDHHHGCIIHDHHFA"));
		r.setReadName("four_way_split_read");
		
		createBAM(input, header, r);
		
		srr.createSupplementaryAlignments(aligner, input, output, output, true, 100);
		List<SAMRecord> list = getRecords(output);
		assertEquals(4, list.size());
	}
	@Test
	@Category(ExternalAlignerTests.class)
	public void realign_full_read_should_not_repeat_primary_alignment() throws IOException {
		ExternalProcessStreamingAligner aligner = new ExternalProcessStreamingAligner(SamReaderFactory.makeDefault(), ExternalAlignerTests.COMMAND_LINE, ExternalAlignerTests.REFERENCE, 4, new IndexedFastaSequenceFile(ExternalAlignerTests.REFERENCE).getSequenceDictionary());
		BufferedReferenceSequenceFile lookup = new BufferedReferenceSequenceFile(ReferenceSequenceFileFactory.getReferenceSequenceFile(ExternalAlignerTests.REFERENCE));
		ProcessingContext pc = new ProcessingContext(new FileSystemContext(testFolder.getRoot(), 500000), ExternalAlignerTests.REFERENCE, lookup, Lists.newArrayList(), getConfig(testFolder.getRoot()));
		SplitReadRealigner srr = new SplitReadRealigner(pc);
		srr.setRealignEntireRecord(true);
		
		SAMFileHeader header = new SAMFileHeader();
		header.setSequenceDictionary(lookup.getSequenceDictionary());
		header.setSortOrder(SortOrder.coordinate);
		
		SAMRecord r = new SAMRecord(header);
		r.setReferenceIndex(0);
		r.setAlignmentStart(1399998);
		r.setCigarString("301S103M");
		r.setReadBases(B("GGATATATAGGGATAGAAGCTTGAATAGTCTGGACATATATTTGTATTGAAATACAAATGTAAGATTTCAGTTAATCAATTTAAACATTTTTATTTTCAAGGGCTTCCAGCGTCCACTTCCTACGGCAAGCAGGAGGAGACAAGCGCCACCCTGCGCTCGCGGAGCCGACCCCGGCTCTCCCCTCCCGTGGCCGCAGGGGTCTGACAGAAAGGGGTCACTAATCTACTTGGCCTTTTGAGGACTGATCCTTAAGAATAATTTTTTTTTTTTTATGATCTTGAAGGCTGAGAAGTATTAGAGTAGGTTTTTTTCTCCTTCATAAGGCCAGATTCTTCTTTCTGTCACAGATTTCAAGTCCCCGCCTCAGCAGCCTTTCACTGTCAGTTCTTTCTCACGTGACCCT"));
		r.setBaseQualities(B("?????BBBB@DEDDDDGGGGGEIEHIHEFHIHIIEHHIEIIIIIIEHII?HHFHHHHDIHIHEHHFIIBCHI=GHIH@HFCEIGIHIDHHHGCIIHDHHFA?????BBBB@DEDDDDGGGGGEIEHIHEFHIHIIEHHIEIIIIIIEHII?HHFHHHHDIHIHEHHFIIBCHI=GHIH@HFCEIGIHIDHHHGCIIHDHHFA?????BBBB@DEDDDDGGGGGEIEHIHEFHIHIIEHHIEIIIIIIEHII?HHFHHHHDIHIHEHHFIIBCHI=GHIH@HFCEIGIHIDHHHGCIIHDHHFA?????BBBB@DEDDDDGGGGGEIEHIHEFHIHIIEHHIEIIIIIIEHII?HHFHHHHDIHIHEHHFIIBCHI=GHIH@HFCEIGIHIDHHHGCIIHDHHFA"));
		r.setReadName("four_way_split_read");
		
		createBAM(input, header, r);
		
		srr.createSupplementaryAlignments(aligner, input, output, output, true, 100);
		List<SAMRecord> list = getRecords(output);
		assertEquals(4, list.size());
	}
	
	@Test
	@Category(ExternalAlignerTests.class)
	public void realign_entire_read_should_drop_existing_supplementary_alignments() throws IOException {
		ExternalProcessStreamingAligner aligner = new ExternalProcessStreamingAligner(SamReaderFactory.makeDefault(), ExternalAlignerTests.COMMAND_LINE, ExternalAlignerTests.REFERENCE, 4, new IndexedFastaSequenceFile(ExternalAlignerTests.REFERENCE).getSequenceDictionary());
		BufferedReferenceSequenceFile lookup = new BufferedReferenceSequenceFile(ReferenceSequenceFileFactory.getReferenceSequenceFile(ExternalAlignerTests.REFERENCE));
		ProcessingContext pc = new ProcessingContext(new FileSystemContext(testFolder.getRoot(), 500000), ExternalAlignerTests.REFERENCE, lookup, Lists.newArrayList(), getConfig(testFolder.getRoot()));
		SplitReadRealigner srr = new SplitReadRealigner(pc);
		srr.setRealignEntireRecord(true);
		SAMFileHeader header = new SAMFileHeader();
		header.setSequenceDictionary(lookup.getSequenceDictionary());
		header.setSortOrder(SortOrder.coordinate);
		
		SAMRecord r = new SAMRecord(header);
		r.setReferenceIndex(0);
		r.setAlignmentStart(1399998);
		r.setCigarString("301S103M");
		r.setReadBases(B("GGATATATAGGGATAGAAGCTTGAATAGTCTGGACATATATTTGTATTGAAATACAAATGTAAGATTTCAGTTAATCAATTTAAACATTTTTATTTTCAAGGGCTTCCAGCGTCCACTTCCTACGGCAAGCAGGAGGAGACAAGCGCCACCCTGCGCTCGCGGAGCCGACCCCGGCTCTCCCCTCCCGTGGCCGCAGGGGTCTGACAGAAAGGGGTCACTAATCTACTTGGCCTTTTGAGGACTGATCCTTAAGAATAATTTTTTTTTTTTTATGATCTTGAAGGCTGAGAAGTATTAGAGTAGGTTTTTTTCTCCTTCATAAGGCCAGATTCTTCTTTCTGTCACAGATTTCAAGTCCCCGCCTCAGCAGCCTTTCACTGTCAGTTCTTTCTCACGTGACCCT"));
		r.setBaseQualities(B("?????BBBB@DEDDDDGGGGGEIEHIHEFHIHIIEHHIEIIIIIIEHII?HHFHHHHDIHIHEHHFIIBCHI=GHIH@HFCEIGIHIDHHHGCIIHDHHFA?????BBBB@DEDDDDGGGGGEIEHIHEFHIHIIEHHIEIIIIIIEHII?HHFHHHHDIHIHEHHFIIBCHI=GHIH@HFCEIGIHIDHHHGCIIHDHHFA?????BBBB@DEDDDDGGGGGEIEHIHEFHIHIIEHHIEIIIIIIEHII?HHFHHHHDIHIHEHHFIIBCHI=GHIH@HFCEIGIHIDHHHGCIIHDHHFA?????BBBB@DEDDDDGGGGGEIEHIHEFHIHIIEHHIEIIIIIIEHII?HHFHHHHDIHIHEHHFIIBCHI=GHIH@HFCEIGIHIDHHHGCIIHDHHFA"));
		r.setReadName("supp_alignment");
		r.setSupplementaryAlignmentFlag(true);
		
		createBAM(input, header, r);
		srr.createSupplementaryAlignments(aligner, input, output, output, true, 100);
		List<SAMRecord> list = getRecords(output);
		assertEquals(0, list.size());	
	}
	@Test
	@Category(ExternalAlignerTests.class)
	public void realign_split_reads_should_drop_existing_supplementary_alignments() throws IOException {
		ExternalProcessStreamingAligner aligner = new ExternalProcessStreamingAligner(SamReaderFactory.makeDefault(), ExternalAlignerTests.COMMAND_LINE, ExternalAlignerTests.REFERENCE, 4, new IndexedFastaSequenceFile(ExternalAlignerTests.REFERENCE).getSequenceDictionary());
		BufferedReferenceSequenceFile lookup = new BufferedReferenceSequenceFile(ReferenceSequenceFileFactory.getReferenceSequenceFile(ExternalAlignerTests.REFERENCE));
		ProcessingContext pc = new ProcessingContext(new FileSystemContext(testFolder.getRoot(), 500000), ExternalAlignerTests.REFERENCE, lookup, Lists.newArrayList(), getConfig(testFolder.getRoot()));
		SplitReadRealigner srr = new SplitReadRealigner(pc);
		srr.setRealignExistingSplitReads(true);
		SAMFileHeader header = new SAMFileHeader();
		header.setSequenceDictionary(lookup.getSequenceDictionary());
		header.setSortOrder(SortOrder.coordinate);
		
		SAMRecord r = new SAMRecord(header);
		r.setReferenceIndex(0);
		r.setAlignmentStart(1399998);
		r.setCigarString("301S103M");
		r.setReadBases(B("GGATATATAGGGATAGAAGCTTGAATAGTCTGGACATATATTTGTATTGAAATACAAATGTAAGATTTCAGTTAATCAATTTAAACATTTTTATTTTCAAGGGCTTCCAGCGTCCACTTCCTACGGCAAGCAGGAGGAGACAAGCGCCACCCTGCGCTCGCGGAGCCGACCCCGGCTCTCCCCTCCCGTGGCCGCAGGGGTCTGACAGAAAGGGGTCACTAATCTACTTGGCCTTTTGAGGACTGATCCTTAAGAATAATTTTTTTTTTTTTATGATCTTGAAGGCTGAGAAGTATTAGAGTAGGTTTTTTTCTCCTTCATAAGGCCAGATTCTTCTTTCTGTCACAGATTTCAAGTCCCCGCCTCAGCAGCCTTTCACTGTCAGTTCTTTCTCACGTGACCCT"));
		r.setBaseQualities(B("?????BBBB@DEDDDDGGGGGEIEHIHEFHIHIIEHHIEIIIIIIEHII?HHFHHHHDIHIHEHHFIIBCHI=GHIH@HFCEIGIHIDHHHGCIIHDHHFA?????BBBB@DEDDDDGGGGGEIEHIHEFHIHIIEHHIEIIIIIIEHII?HHFHHHHDIHIHEHHFIIBCHI=GHIH@HFCEIGIHIDHHHGCIIHDHHFA?????BBBB@DEDDDDGGGGGEIEHIHEFHIHIIEHHIEIIIIIIEHII?HHFHHHHDIHIHEHHFIIBCHI=GHIH@HFCEIGIHIDHHHGCIIHDHHFA?????BBBB@DEDDDDGGGGGEIEHIHEFHIHIIEHHIEIIIIIIEHII?HHFHHHHDIHIHEHHFIIBCHI=GHIH@HFCEIGIHIDHHHGCIIHDHHFA"));
		r.setReadName("supp_alignment");
		r.setSupplementaryAlignmentFlag(true);
		
		createBAM(input, header, r);
		srr.createSupplementaryAlignments(aligner, input, output, output, true, 100);
		List<SAMRecord> list = getRecords(output);
		assertEquals(0, list.size());
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
		SplitReadRealigner srr = new SplitReadRealigner(pc);
		srr.setRealignEntireRecord(true);
		srr.setReference(ref);
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
		srr.createSupplementaryAlignments(aligner, input, output, output, true, 1);
		List<SAMRecord> list = getRecords(output);
		assertEquals(2, list.size());
		Assert.assertEquals("282S618M", list.get(0).getCigarString());
	}
	//@Test // presort=false was the issue!
	public void mergeSupplementaryAlignment_debug_memory_leak() throws IOException {
		File projectDir = new File("W:/projects/BouilletProject1/working");
		File asmDir = new File(projectDir, "bouillet.asm.bam.gridss.working");
		ProcessingContext pc = new ProcessingContext(
				new FileSystemContext(projectDir, 500000),
				new File("W:/projects/BouilletProject1/ref/mm10.fa"),
				null,
				Lists.newArrayList(),
				new GridssConfiguration(getDefaultConfig(), projectDir));
		SplitReadRealigner srr = new SplitReadRealigner(pc);
		List<File> realignments = ImmutableList.of(
				new File(asmDir, "bouillet.asm.bam.realign.0.bam"),
				new File(asmDir, "bouillet.asm.bam.realign.1.bam"),
				new File(asmDir, "bouillet.asm.bam.realign.2.bam"),
				new File(asmDir, "bouillet.asm.bam.realign.3.bam"),
				new File(asmDir, "bouillet.asm.bam.realign.4.bam"),
				new File(asmDir, "bouillet.asm.bam.realign.5.bam"),
				new File(asmDir, "bouillet.asm.bam.realign.6.bam"),
				new File(asmDir, "bouillet.asm.bam.realign.7.bam")
		);
		srr.mergeSupplementaryAlignment(
				new File("W:/projects/BouilletProject1/bouillet.asm.bam"),
				realignments,
				new File(projectDir, "debug_memory_leak_order.bam"),
				new File(projectDir, "debug_memory_leak_unordered.bam"),
				true);
	}
}
