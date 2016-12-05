package gridss;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.IntermediateFilesTest;
import au.edu.wehi.idsv.NonReferenceReadPair;
import au.edu.wehi.idsv.SoftClipEvidence;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;

public class CallVariantsTest extends IntermediateFilesTest {
	@Test
	public void test_sr_rp_assembly() throws IOException {
		List<DirectedEvidence> in = new ArrayList<DirectedEvidence>();
		in.add(SCE(FWD, withSequence("AACCGGTTCTA", Read(0, 15, "5M6S"))));
		in.add(SCE(FWD, withSequence("AACCGGTTCTA", Read(0, 15, "6M5S"))));
		in.add(NRRP(withSequence("AACCGGTTCTA", DP(0, 1, "11M", true, 1, 100, "11M", false))));
		in.add(NRRP(withSequence("AACCGGTTCTA", DP(0, 2, "11M", true, 1, 100, "11M", false))));
		
		List<SAMRecord> insam = in.stream().flatMap(de -> {
			if (de instanceof SoftClipEvidence) return ImmutableList.<SAMRecord>of(((SoftClipEvidence)de).getSAMRecord()).stream();
			if (de instanceof NonReferenceReadPair) return ImmutableList.<SAMRecord>of(
					((NonReferenceReadPair)de).getNonReferenceRead(),
					((NonReferenceReadPair)de).getLocalledMappedRead()).stream();
			throw new RuntimeException("NYI");
		}).collect(Collectors.toList());
		insam.addAll(Lists.newArrayList(RP(0, 1, 100, 10)));
		createInput(insam);
		File propfile = testFolder.newFile("custom.properties");
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(propfile))) {
			writer.write("assembly.k=5\n");
			writer.write("softclip.minLength=1\n");
			writer.write("minAnchorShannonEntropy=0\n");
			writer.write("softclip.minAnchorIdentity=0\n");
			writer.write("softclip.minAverageQual=0\n");
			writer.write("realignment.aligner=\n");
		}
		File assembly = new File(testFolder.getRoot(), "assembly.bam");
		createBAM(assembly, SortOrder.coordinate);
		String[] args = new String[] {
				"INPUT=" + input.toString(),
				"INPUT_CATEGORY=5",
				"ASSEMBLY=" + assembly.toString(),
				"REFERENCE_SEQUENCE=" + reference.toString(),
				"OUTPUT=" + output.toString(),
				"TMP_DIR=" + super.testFolder.getRoot().toString(),
				"WORKING_DIR=" + super.testFolder.getRoot().toString(),
				"CONFIGURATION_FILE=" + propfile.toString(),
				"INPUT_MIN_FRAGMENT_SIZE=10",
				"INPUT_MAX_FRAGMENT_SIZE=100",
		};
		assertEquals(0, new CallVariants().instanceMain(args));
		List<SAMRecord> breakendAssemblies = getRecords(new File(output.getAbsoluteFile() + ".gridss.working/" + output.getName() + ".assembly.bam"));
		assertEquals(1, breakendAssemblies.size());
		assembly.delete();
	}
	@Test
	public void should_handle_unpaired_libraries() throws IOException {
		List<DirectedEvidence> in = new ArrayList<DirectedEvidence>();
		in.add(SCE(FWD, withSequence("AACCGGTTCTA", Read(0, 15, "5M6S"))));
		
		List<SAMRecord> insam = in.stream().flatMap(de -> {
			if (de instanceof SoftClipEvidence) return ImmutableList.<SAMRecord>of(((SoftClipEvidence)de).getSAMRecord()).stream();
			if (de instanceof NonReferenceReadPair) return ImmutableList.<SAMRecord>of(
					((NonReferenceReadPair)de).getNonReferenceRead(),
					((NonReferenceReadPair)de).getLocalledMappedRead()).stream();
			throw new RuntimeException("NYI");
		}).collect(Collectors.toList());
		createInput(insam);
		File assembly = new File(testFolder.getRoot(), "assembly.bam");
		createBAM(assembly, SortOrder.coordinate);
		String[] args = new String[] {
				"INPUT=" + input.toString(),
				"ASSEMBLY=" + assembly.toString(),
				"REFERENCE_SEQUENCE=" + reference.toString(),
				"OUTPUT=" + output.toString(),
				"TMP_DIR=" + super.testFolder.getRoot().toString(),
				"WORKING_DIR=" + super.testFolder.getRoot().toString()
		};
		assertEquals(0, new CallVariants().instanceMain(args));
		assertTrue(output.exists());
		assembly.delete();
	}
}
