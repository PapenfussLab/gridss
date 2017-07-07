package gridss;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;

import com.google.common.collect.Lists;

import au.edu.wehi.idsv.IntermediateFilesTest;
import au.edu.wehi.idsv.bed.BedpeIterator;
import au.edu.wehi.idsv.bed.BedpeRecord;

public class AnnotateInexactHomologyBedpeTest extends IntermediateFilesTest {
	@Test
	public void should_set_score_to_homology_length() throws FileNotFoundException {
		AnnotateInexactHomologyBedpe cmd = new AnnotateInexactHomologyBedpe();
		Assert.assertEquals(0, cmd.instanceMain(new String[] {
				"R=" + SMALL_FA_FILE.getAbsolutePath(),
				"I=" + new File("src/test/resources/simple.bedpe"),
				"O=" + output,
				"M=0",
				"D=5",
		}));
		List<BedpeRecord> result = Lists.newArrayList(new BedpeIterator(output, SMALL_FA));
		Assert.assertEquals(3, result.size());
		Assert.assertEquals(10, result.get(1).score, 0);
	}
	@Test
	public void should_include_untemplated_sequence_in_homology_calculation() throws FileNotFoundException {
		AnnotateInexactHomologyBedpe cmd = new AnnotateInexactHomologyBedpe();
		Assert.assertEquals(0, cmd.instanceMain(new String[] {
				"R=" + SMALL_FA_FILE.getAbsolutePath(),
				"I=" + new File("src/test/resources/simple.bedpe"),
				"O=" + output,
				"M=10",
				"D=10",
				"UC=11",
		}));
		List<BedpeRecord> result = Lists.newArrayList(new BedpeIterator(output, SMALL_FA));
		Assert.assertEquals(3, result.size()); 
		Assert.assertEquals(5, result.get(2).score, 0);
	}
}
