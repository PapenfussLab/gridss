package au.edu.wehi.idsv;

import static org.junit.Assert.assertEquals;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;

import java.util.Collections;
import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.vcf.VcfAttributes.Subset;



public class SequentialCoverageAnnotatorTest extends TestHelper {
	public VariantContextDirectedEvidence go(List<DirectedEvidence> evidence, List<SAMRecord> ref, VariantContextDirectedEvidence toAnnotate) {
		Collections.sort(ref, new SAMRecordCoordinateComparator());
		return new SequentialCoverageAnnotator(getContext(), new SequentialReferenceCoverageLookup(ref.iterator(), 1024), null).annotate(toAnnotate);
	}
	@Test
	public void should_add_reference_counts() {
		VariantContextDirectedEvidence result = go(L(
				(DirectedEvidence)SoftClipEvidence.create(getContext(), SES(), FWD, withSequence("TTTTTTTT", Read(0, 1, "3S1M3S"))[0], withSequence("TTT", Read(0, 12, "3M"))[0])),
			L(
				RP(0, 1, 100, 5),
				RP(0, 2, 100, 5),
				RP(0, 9, 100, 10)),
			(VariantContextDirectedEvidence)minimalBreakend()
				.breakpoint(new BreakpointSummary(0, FWD, 10, 10, 1, BWD, 100, 100), "")
				.make());
		assertEquals(1, result.getReferenceReadCount(Subset.NORMAL));
		assertEquals(2, result.getReferenceReadPairCount(Subset.NORMAL));
	}
}
