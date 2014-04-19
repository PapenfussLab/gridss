package au.edu.wehi.socrates;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;

import org.junit.Test;


public class FileNamingConventionTest {
	@Test
	public void getMetrics_should_use_suffix() throws IOException {
		assertEquals(new File("test.bam.socrates.insertsize.txt").getAbsolutePath(), FileNamingConvention.GetMetrics(new File("test.bam")).getAbsolutePath());
	}
	@Test
	public void getMetrics_should_return_common_metrics_for_all_associated_files() throws IOException {
		assertEquals(new File("test.bam.socrates.insertsize.txt").getAbsolutePath(), FileNamingConvention.GetMetrics(new File("test.bam.socrates.sv.bam")).getAbsolutePath());
	}
}
