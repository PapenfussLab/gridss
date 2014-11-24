package au.edu.wehi.idsv.visualisation;

import htsjdk.samtools.util.Log;

import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.File;
import java.io.FileWriter;
import java.util.PriorityQueue;

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
	private PriorityQueue<SubgraphAssemblyAlgorithmTracker> sortBuffer = new PriorityQueue<SubgraphAssemblyAlgorithmTracker>(32, SubgraphAssemblyAlgorithmTracker.ByGenomicPosition);
	private int windowSize;
	public SubgraphAssemblyAlgorithmTrackerBEDWriter(int maxSubgraphSize, File bedFile) {
		this.windowSize = 2 * maxSubgraphSize + 2;
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
		if (tracker == null) return;
		if (tracker.getReferenceIndex() < 0) return;
		if (tracker.getStartAnchorPosition() < 0) return;
		if (tracker.getEndAnchorPosition() < 0) return;
		sortBuffer.add(tracker);
		flushPosition(tracker.getReferenceIndex(), tracker.getStartAnchorPosition() - windowSize - 1);
	}
	private void flushPosition(int referenceIndex, long pos) {
		while (!sortBuffer.isEmpty() && (sortBuffer.peek().getReferenceIndex() < referenceIndex || sortBuffer.peek().getStartAnchorPosition() <= pos)) {
			SubgraphAssemblyAlgorithmTracker tracker = sortBuffer.poll();
			try {
				if (bw != null) {
					String result = tracker.toBed();
					if (!StringUtils.isEmpty(result)) {
						bw.append(result);
					}
				}
			} catch (Exception e) {
				log.warn("Unable to write de Bruijn algorithm tracking data", e);
				close(false);
			}
		}
	}
	public void close(boolean flush) {
		if (flush) {
			flushPosition(Integer.MAX_VALUE, Long.MAX_VALUE);
		}
		try {
			if (bw != null) bw.close();
			bw = null;
			if (fw != null) fw.close();
			fw = null;
		} catch (Exception e) {
			log.warn("Error closing de Bruijn algorithm tracker stream", e);
		}
	}
	@Override
	public void close() {
		close(true);
	}
}
