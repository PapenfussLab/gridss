package performancetesting;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;

import java.io.File;

import org.junit.Ignore;
import org.junit.Test;

import au.edu.wehi.idsv.TestHelper;

public class BamReadTest extends TestHelper {
	@Ignore // TODO: performance testing group
	@Test
	public void read() {
		//String file = "W:/i/archive/data.na12878/00000000000000000000000000000000.sc.bam";
		String file = "W:/projects/DREAM/download/PCSI_0072_Si_R_PE_762_WG_NoIndex_6_130905_SN801_0125_AC2D9BACXX.bam";
		CloseableIterator<SAMRecord> it = getContext().getSamReaderIterator(new File(file));
		while (it.hasNext()) it.next();
		it.close();
	}
}
