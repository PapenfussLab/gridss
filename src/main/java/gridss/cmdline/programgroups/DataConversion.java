package gridss.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;

public class DataConversion implements CommandLineProgramGroup {
	@Override
    public String getName() { return "Data Conversion"; }
    @Override
    public String getDescription() { return "Tools for changing data representation formats."; }
}
