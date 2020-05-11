package au.edu.wehi.idsv.model;

import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;
import gridss.analysis.CigarDetailMetrics;
import gridss.analysis.IdsvMetrics;
import gridss.analysis.InsertSizeDistribution;
import gridss.analysis.MapqMetrics;
import org.junit.Test;
import picard.analysis.InsertSizeMetrics;

import java.util.ArrayList;

import static org.junit.Assert.assertEquals;

public class EmpiricalReferenceLikelihoodModelTest extends TestHelper {
	private static final EmpiricalReferenceLikelihoodModel model = new EmpiricalReferenceLikelihoodModel();
	@Test
	public void should_score_read_pairs_correctly() {
		IdsvSamFileMetrics metrics = new IdsvSamFileMetrics(
				new InsertSizeMetrics() {{
					MEAN_INSERT_SIZE = 2;
					MEDIAN_INSERT_SIZE = 2;
					MIN_INSERT_SIZE = 1;
					MAX_INSERT_SIZE = 3;
					MEDIAN_ABSOLUTE_DEVIATION = 1;
				}},
				new IdsvMetrics() {{
					MAX_READ_LENGTH = 100;
					MAX_PROPER_PAIR_FRAGMENT_LENGTH = 2;
					MIN_PROPER_PAIR_FRAGMENT_LENGTH = 2;
					READ_PAIRS = 13;
					READ_PAIRS_ONE_MAPPED = 2;
					READ_PAIRS_ZERO_MAPPED = 1;
					READ_PAIRS_BOTH_MAPPED = READ_PAIRS - READ_PAIRS_ONE_MAPPED - READ_PAIRS_ZERO_MAPPED;
					READS = 2 * READ_PAIRS;
					MAPPED_READS = READS - READ_PAIRS_ONE_MAPPED - 2*READ_PAIRS_ZERO_MAPPED;
				}},
				new MapqMetrics() {{
				}},
				new InsertSizeDistribution(
					new int[] { 1, 2, 3,},
					new double[] { 1, 8, 1, }),
				new ArrayList<CigarDetailMetrics>());
		assertEquals(10, model.scoreReadPair(metrics, 3, 1000, 1000), 0.001);
	}
}
