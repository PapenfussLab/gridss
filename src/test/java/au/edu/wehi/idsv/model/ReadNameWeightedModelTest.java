package au.edu.wehi.idsv.model;

import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.TestHelper;
import org.junit.Test;

import static org.junit.Assert.assertEquals;


public class ReadNameWeightedModelTest extends TestHelper {
	private static final ReadNameWeightedModel model = new ReadNameWeightedModel("regex_r_(?<weight>[0-9]+)_");
	@Test
	public void should_fall_back_to_read_count() {
		DirectedEvidence e = SCE(FWD, Read(0, 1, "50M50S"));
		assertEquals(1, model.scoreIndel(null, e, null, -1, 5), 0);
		assertEquals(1, model.scoreSoftClip(null,e, -1, 5), 0);
		assertEquals(1, model.scoreUnmappedMate(null,e, 5), 0);
		assertEquals(1, model.scoreSplitRead(null,e,  0, 10, 10), 0.000001);
		assertEquals(1, model.scoreReadPair(null,e,  0, 10, 10), 0.000001);
	}

	/**
	 * Assembly quality is a pass-though of the underlying scoring model
	 */
	@Test
	public void asm_should_fall_back_to_asm_qual() {
		DirectedEvidence e = SCE(FWD, Read(0, 1, "50M50S"));
		assertEquals(10, model.scoreAssembly(e, 1, 2, 4, 8, 16, 32), 0.000001);
		assertEquals(10, model.scoreBreakendAssembly(e, 1, 2, 4, 8, 16), 0.000001);
	}
	@Test
	public void use_evidence_if_regex_matches() {
		DirectedEvidence e = SCE(FWD, withReadName("regex_r_69_", Read(0, 1, "50M50S"))[0]);
		assertEquals(69, model.scoreIndel(null, e, null, -1, 5), 0);
		assertEquals(69, model.scoreSoftClip(null,e, -1, 5), 0);
		assertEquals(69, model.scoreUnmappedMate(null,e, 5), 0);
		assertEquals(69, model.scoreSplitRead(null,e,  0, 10, 10), 0.000001);
		assertEquals(69, model.scoreReadPair(null,e,  0, 10, 10), 0.000001);
		assertEquals(69, model.scoreAssembly(e, 1, 2, 4, 8, 16, 32), 0.000001);
		assertEquals(69, model.scoreBreakendAssembly(e, 1, 2, 4, 8, 16), 0.000001);
	}
}
