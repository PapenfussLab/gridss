package au.edu.wehi.idsv.util;

import au.edu.wehi.idsv.TestHelper;
import org.junit.Assert;
import org.junit.Test;

public class FilenameUtilTest extends TestHelper {
    @Test
    public void shouldWorkForVirusBreakendTelemetry() {
        Assert.assertEquals("krakentaxid129951NC_001405", FilenameUtil.stripInvalidFilenameCharacters("kraken:taxid|129951|NC_001405"));
    }
}