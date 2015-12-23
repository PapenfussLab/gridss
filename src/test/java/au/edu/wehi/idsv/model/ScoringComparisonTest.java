package au.edu.wehi.idsv.model;

import java.util.List;

import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;
import au.edu.wehi.idsv.metrics.CigarDetailMetrics;
import au.edu.wehi.idsv.metrics.CigarSizeDistributionTest;
import au.edu.wehi.idsv.metrics.IdsvSamFileMetrics;

public class ScoringComparisonTest extends TestHelper {
	private static int[] mapq = new int[] { 0, 5, 10, 15, 20, 25, 30, 35, 40 };
	private void printSoftClipComparison(VariantScoringModel model1, VariantScoringModel model2) {
		List<CigarDetailMetrics> sc = CigarSizeDistributionTest.data_778();
		IdsvSamFileMetrics metrics = new IdsvSamFileMetrics(null, null, null, sc);
		System.out.println("Split Read:");
		System.out.println("mapq1,mapq2,sclength,score1,score2");
		for (int i = 0; i < mapq.length; i++) {
			int mapq1 = mapq[i];
			for (int j = i; j < mapq.length; j++) {
				int mapq2 = mapq[j];
				for (int scLength = 0; scLength < sc.size(); scLength += 5) {
					System.out.println(String.format("%d,%d,%d,%f,%f", mapq1, mapq2, scLength,
							model1.scoreSplitRead(metrics, scLength, mapq1, mapq2),
							model2.scoreSplitRead(metrics, scLength, mapq1, mapq2)));
				}
			}
		}
		System.out.println("SoftClip:");
		System.out.println("mapq1,sclength,score1,score2");
		for (int i = 0; i < mapq.length; i++) {
			int mapq1 = mapq[i];
			for (int scLength = 0; scLength < sc.size(); scLength += 5) {
				System.out.println(String.format("%d,%d,%f,%f", mapq1, scLength,
						model1.scoreSoftClip(metrics, scLength, mapq1),
						model2.scoreSoftClip(metrics, scLength, mapq1)));
			}
		}
	}
	@Test
	public void SoftClip_Fast_Empirical() {
		printSoftClipComparison(new FastEmpiricalReferenceLikelihoodModel(), new EmpiricalLlrModel());
	}
	@Test
	public void Fast_Error() {
		printSoftClipComparison(new FastEmpiricalReferenceLikelihoodModel(), new EmpiricalReferenceLikelihoodModel());
	}
}
