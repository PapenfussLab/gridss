package au.edu.wehi.idsv.debruijn.positional;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.DirectedEvidenceOrder;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.TestHelper;

import com.google.common.collect.Lists;


public class PositionalAssemblerTest extends TestHelper {
	@Test
	public void should_assemble_simple_input() {
		ProcessingContext pc = getContext();
		AssemblyEvidenceSource aes = AES(pc);
		pc.getAssemblyParameters().k = 4;
		List<DirectedEvidence> input = new ArrayList<DirectedEvidence>();
		input.add(SCE(BWD, Read(0, 10, "5S5M")));
		input.add(SCE(FWD, Read(0, 10, "5M5S")));
		input.add(SCE(BWD, Read(0, 100, "5S5M")));
		input.add(SCE(FWD, Read(0, 100, "5M5S")));
		input.sort(DirectedEvidenceOrder.ByStartEnd);
		ArrayList<SAMRecordAssemblyEvidence> r = Lists.newArrayList(new PositionalAssembler(pc, aes, input.iterator()));
		assertEquals(4, r.size());
	}
}
