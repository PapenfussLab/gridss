package au.edu.wehi.idsv;

import java.io.IOException;

import org.junit.Assert;
import org.junit.Test;

import com.google.common.io.Files;

public class PrecomputedGcBiasAdjusterTest extends IntermediateFilesTest {
	@Test
	public void should_require_0_to_100_multiplers_to_be_set() throws IOException {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i <= 100; i++) {
			sb.append(i);
			sb.append('\t');
			sb.append(i + 0.5);
			sb.append('\n');
		}
		Files.write(B(sb.toString()), input);
		PrecomputedGcBiasAdjuster gcba = new PrecomputedGcBiasAdjuster(input);
		for (int i = 0; i <= 100; i++) {
			Assert.assertEquals(i + 0.5, gcba.adjustmentMultiplier(i), 0);
		}
	}
	@Test(expected=IllegalArgumentException.class)
	public void should_require_all_values_set() throws IOException {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < 100; i++) { // missing adjustment for 100% GC 
			sb.append(i);
			sb.append('\t');
			sb.append(i + 0.5);
			sb.append('\n');
		}
		Files.write(B(sb.toString()), input);
		new PrecomputedGcBiasAdjuster(input);
	}
}
