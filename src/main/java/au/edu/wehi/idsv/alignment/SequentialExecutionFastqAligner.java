package au.edu.wehi.idsv.alignment;

import java.io.File;
import java.io.IOException;

/**
 * Fastq aligner in which only one alignment can be performed at any given time.
 * As aligners typically use all available cores, throughput can be improved by
 * ensuring only one instance is running at any given time.
 */
public class SequentialExecutionFastqAligner implements FastqAligner {
	private static Object lock = new Object();
	private final FastqAligner aligner;
	public SequentialExecutionFastqAligner(FastqAligner aligner) {
		this.aligner = aligner;
	}
	public void align(File fastq, File output, File reference, int threads) throws IOException {
		synchronized(lock) {
			aligner.align(fastq, output, reference, threads);
		}
	}
}
