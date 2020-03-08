package au.edu.wehi.idsv;

import org.junit.Assert;
import org.junit.Test;

public class ReadPairConcordanceCalculatorTest {
	@Test
	public void create_should_create_calculator() {
		Assert.assertNotNull(ReadPairConcordanceCalculator.create(0, 100, null, null, null));
	}
}
