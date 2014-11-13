package au.edu.wehi.idsv.visualisation;

import htsjdk.samtools.util.Log;

import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.File;
import java.io.FileWriter;

import org.apache.commons.lang3.StringUtils;

/**
 * Writes algorithm tracking information to a BED file
 * @author Daniel Cameron
 *
 */
public class SubgraphAssemblyAlgorithmTrackerBEDWriter implements Closeable {
	private static final Log log = Log.getInstance(SubgraphAssemblyAlgorithmTrackerBEDWriter.class);
	private FileWriter fw = null;
	private BufferedWriter bw = null;
	public SubgraphAssemblyAlgorithmTrackerBEDWriter(File bedFile) {
		try {
			fw = new FileWriter(bedFile.getAbsoluteFile());
			bw = new BufferedWriter(fw);
			//bw.append("##");
			bw.append("track gffTags=on\n");
		} catch (Exception e) {
			log.warn(e, "Unable to track progress to " + bedFile.toString());
			close();
		}
	}
	public void write(SubgraphAssemblyAlgorithmTracker tracker) {
		try {
			if (bw != null) {
				String result = tracker.toBed();
				if (!StringUtils.isEmpty(result)) {
					bw.append(result);
				}
			}
		} catch (Exception e) {
			log.debug("Unable to write de Bruijn algorithm tracking data", e);
		}
	}
	@Override
	public void close() {
		try {
			if (bw != null) bw.close();
			bw = null;
			if (fw != null) fw.close();
			fw = null;
		} catch (Exception e) {
			log.debug("Error closing de Bruijn algorithm tracker stream", e);
		}
	}
}
