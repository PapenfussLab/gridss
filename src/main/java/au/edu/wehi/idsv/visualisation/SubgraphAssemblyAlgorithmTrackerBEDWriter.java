package au.edu.wehi.idsv.visualisation;

import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.File;
import java.io.FileWriter;
import java.util.PriorityQueue;

import org.apache.commons.lang3.StringUtils;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;

import au.edu.wehi.idsv.graph.PathNode;
import htsjdk.samtools.util.Log;

/**
 * Writes algorithm tracking information to a BED file
 * @author Daniel Cameron
 *
 */
public class SubgraphAssemblyAlgorithmTrackerBEDWriter<T, PN extends PathNode<T>> implements Closeable {
	private static final Log log = Log.getInstance(SubgraphAssemblyAlgorithmTrackerBEDWriter.class);
	public Ordering<SubgraphAssemblyAlgorithmTracker<T, PN>> ByGenomicPosition = new Ordering<SubgraphAssemblyAlgorithmTracker<T, PN>>() {
		@Override
		public int compare(SubgraphAssemblyAlgorithmTracker<T, PN> left, SubgraphAssemblyAlgorithmTracker<T, PN> right) {
			return ComparisonChain.start()
			        .compare(left.getReferenceIndex(), right.getReferenceIndex())
			        .compare(left.getStartAnchorPosition(), right.getStartAnchorPosition())
			        .compare(left.getEndAnchorPosition(), right.getEndAnchorPosition())
			        .result();
		}
	};
	private FileWriter fw = null;
	private BufferedWriter bw = null;
	private PriorityQueue<SubgraphAssemblyAlgorithmTracker<T, PN>> sortBuffer = new PriorityQueue<SubgraphAssemblyAlgorithmTracker<T, PN>>(32, ByGenomicPosition);
	private int windowSize;
	public SubgraphAssemblyAlgorithmTrackerBEDWriter(int maxSubgraphSize, File bedFile) {
		this.windowSize = 16 * maxSubgraphSize + 1024;
		try {
			bedFile.getParentFile().mkdirs();
			fw = new FileWriter(bedFile.getAbsoluteFile());
			bw = new BufferedWriter(fw);
			//bw.append("##");
			bw.append("track gffTags=on useScore=1 \n");
		} catch (Exception e) {
			log.warn(e, "Unable to track progress to " + bedFile.toString());
			close();
		}
	}
	public void write(SubgraphAssemblyAlgorithmTracker<T, PN> tracker) {
		if (tracker == null) return;
		if (tracker.getReferenceIndex() < 0) return;
		if (tracker.getStartAnchorPosition() < 0) return;
		if (tracker.getEndAnchorPosition() < 0) return;
		sortBuffer.add(tracker);
		flushPosition(tracker.getReferenceIndex(), tracker.getStartAnchorPosition() - windowSize - 1);
	}
	private void flushPosition(int referenceIndex, long pos) {
		while (!sortBuffer.isEmpty() && (sortBuffer.peek().getReferenceIndex() < referenceIndex || sortBuffer.peek().getStartAnchorPosition() <= pos)) {
			SubgraphAssemblyAlgorithmTracker<T, PN> tracker = sortBuffer.poll();
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
