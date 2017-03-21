package au.edu.wehi.idsv;

import org.junit.Assert;
import org.junit.Test;

public class ReadPairConcordanceCalculatorTest {
	@Test
	public void create_should_create_calculator() {
		Assert.assertNotNull(ReadPairConcordanceCalculator.create(ReadPairConcordanceMethod.SAM_FLAG, 0, 0, 0, null, null));
	}
}
