package gridss.filter;

import au.edu.wehi.idsv.TestHelper;
import htsjdk.samtools.SAMRecord;
import org.junit.Assert;
import org.junit.Test;

public class SplitReadFilterTest extends TestHelper {
	@Test
	public void should_use_SA_tag() {
		SAMRecord r = Read(0, 1, "10M");
		Assert.assertTrue(new SplitReadFilter().filterOut(r));
		r.setAttribute("SA", "something");
		Assert.assertFalse(new SplitReadFilter().filterOut(r));
	}
}
