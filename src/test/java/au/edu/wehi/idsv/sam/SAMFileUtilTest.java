package au.edu.wehi.idsv.sam;

import au.edu.wehi.idsv.IntermediateFilesTest;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Ordering;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class SAMFileUtilTest extends IntermediateFilesTest {
	@Test
	public void merge_should_merge_in_order() throws IOException {
		File input2 = testFolder.newFile("input2.bam");
		input2.delete();
		createBAM(input, SortOrder.coordinate,
				Read(0, 1, "1M"),
				Read(0, 3, "1M"),
				Read(0, 5, "1M"));
		createBAM(input2, SortOrder.coordinate,
				Read(0, 2, "1M"),
				Read(0, 4, "1M"),
				Read(0, 6, "1M"));
		SAMFileUtil.merge(ImmutableList.of(input, input2), output);
		List<SAMRecord> list = getRecords(output);
		assertEquals(6, list.size());
		assertTrue(Ordering.from(SortOrder.coordinate.getComparatorInstance()).isOrdered(list));
	}
	@Test
	public void sort_should_match_output_order() throws IOException {
		File output = testFolder.newFile("output.bam");
		createBAM(input, SortOrder.unsorted,
				Read(5, 1, "1M"),
				Read(3, 3, "1M"),
				Read(1, 5, "1M"));
		SAMFileUtil.sort(getFSContext(), input, output, SortOrder.coordinate);
		assertTrue(Ordering.from(SortOrder.coordinate.getComparatorInstance()).isOrdered(getRecords(output)));
	}
	@Test
	public void sort_should_copy_file_for_matching_sort_order() throws IOException {
		File output = testFolder.newFile("output.bam");
		createBAM(input, SortOrder.queryname,
				withReadName("1", Read(5, 1, "1M"))[0],
				withReadName("2", Read(3, 3, "1M"))[0],
				withReadName("2", Read(1, 5, "1M"))[0]);
		SAMFileUtil.sort(getFSContext(), input, output, SortOrder.coordinate);
		assertTrue(Ordering.from(SortOrder.coordinate.getComparatorInstance()).isOrdered(getRecords(output)));
	}
	@Test(expected=IllegalArgumentException.class)
	public void merge_should_require_same_sort_order() throws IOException {
		createBAM(input, SortOrder.queryname,
				withReadName("1", Read(5, 1, "1M"))[0],
				withReadName("2", Read(3, 3, "1M"))[0],
				withReadName("2", Read(1, 5, "1M"))[0]);
		createBAM(output, SortOrder.coordinate,
				withReadName("1", Read(5, 1, "1M"))[0],
				withReadName("2", Read(3, 3, "1M"))[0],
				withReadName("2", Read(1, 5, "1M"))[0]);
		SAMFileUtil.merge(ImmutableList.of(input, output), output);
	}
}
