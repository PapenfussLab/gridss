package gridss;

import au.edu.wehi.idsv.IntermediateFilesTest;
import au.edu.wehi.idsv.bed.BedpeIterator;
import au.edu.wehi.idsv.bed.BedpeRecord;
import com.google.common.collect.Lists;
import gridss.cmdline.CommandLineProgramHelper;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;

public class AnnotateInexactHomologyBedpeTest extends IntermediateFilesTest {
	@Test
	public void should_set_score_to_homology_length() throws FileNotFoundException {
		CommandLineProgramHelper cmd = new CommandLineProgramHelper(new AnnotateInexactHomologyBedpe());
		cmd.addArg("R", SMALL_FA_FILE.getAbsolutePath());
		cmd.addArg("I", new File("src/test/resources/simple.bedpe"));
		cmd.addArg("O", output);
		cmd.addArg("M", 0);
		cmd.addArg("D", 5);
		Assert.assertEquals(0, cmd.run());
		List<BedpeRecord> result = Lists.newArrayList(new BedpeIterator(output, SMALL_FA.getSequenceDictionary()));
		Assert.assertEquals(3, result.size());
		Assert.assertEquals("10", result.get(1).score);
	}
	@Test
	public void should_exclude_untemplated_sequence_from_homology_calculation() throws FileNotFoundException {
		CommandLineProgramHelper cmd = new CommandLineProgramHelper(new AnnotateInexactHomologyBedpe());
		cmd.addArg("R", SMALL_FA_FILE.getAbsolutePath());
		cmd.addArg("I", new File("src/test/resources/simple.bedpe"));
		cmd.addArg("O", output);
		cmd.addArg("M", 10);
		cmd.addArg("D", 10);
		cmd.addArg("UC", 11);
		Assert.assertEquals(0, cmd.run());
		List<BedpeRecord> result = Lists.newArrayList(new BedpeIterator(output, SMALL_FA.getSequenceDictionary()));
		Assert.assertEquals(3, result.size()); 
		Assert.assertEquals("0", result.get(2).score);
	}
}
