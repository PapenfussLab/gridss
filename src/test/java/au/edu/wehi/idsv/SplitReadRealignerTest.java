package au.edu.wehi.idsv;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.junit.Before;
import org.junit.Test;

import com.google.common.collect.Ordering;

import au.edu.wehi.idsv.alignment.AlignerFactory;
import au.edu.wehi.idsv.alignment.SmithWatermanFastqAligner;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;

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
	public void should_not_call_aligner_unless_required() throws IOException {
		createInput(Read(0, 1, "50M"));
		SplitReadRealigner srr = new SplitReadRealigner(getContext(), null);
		srr.createSupplementaryAlignments(input, output);
		assertEquals(1, getRecords(output).size());
	}
	@Test
	public void should_output_sorted_merge() throws IOException {
		SAMRecord r = Read(2, 1, "50S50M");
		r.setReadBases(B(S(RANDOM).substring(100, 150) + S(RANDOM).substring(0, 50)));
		r.setReadName("r");
		
		createBAM(input, SortOrder.coordinate, r);
		SplitReadRealigner srr = new SplitReadRealigner(getContext(), aligner);
		srr.createSupplementaryAlignments(input, output);
		
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
		SplitReadRealigner srr = new SplitReadRealigner(getContext(), aligner);
		srr.createSupplementaryAlignments(input, output);
		
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
		srr = new SplitReadRealigner(getContext(), aligner);
		srr.createSupplementaryAlignments(input, output);
		result = getRecords(output);
		assertEquals(3, result.size());
		assertTrue(Ordering.from(SortOrder.coordinate.getComparatorInstance()).isOrdered(result));
		
		createBAM(input, SortOrder.queryname, r);
		srr = new SplitReadRealigner(getContext(), aligner);
		srr.createSupplementaryAlignments(input, output);
		result = getRecords(output);
		assertEquals(3, result.size());
		assertTrue(Ordering.from(SortOrder.queryname.getComparatorInstance()).isOrdered(result));
		
		createBAM(input, SortOrder.unsorted, r);
		srr = new SplitReadRealigner(getContext(), aligner);
		srr.createSupplementaryAlignments(input, output);
		result = getRecords(output);
		assertEquals(3, result.size());
	}
	@Test
	public void fastq_keys_should_be_unique() throws IOException {
		createBAM(input, SortOrder.coordinate, Read(0, 1, "1S1M1S"));
		new SplitReadRealigner(getContext(), aligner).createSupplementaryAlignments(input, output);
		File fq1 = getContext().getFileSystemContext().getRealignmentFastq(input, 0);
		List<FastqRecord> fq = getFastqRecords(fq1);
		assertEquals(2, fq.size());
		assertNotEquals(fq.get(0).getReadHeader(), fq.get(1).getReadHeader());
	}
}
