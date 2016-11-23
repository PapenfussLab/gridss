package au.edu.wehi.idsv.alignment;

import java.io.File;
import java.io.IOException;

/**
 * Fastq aligner in which only one alignment can be performed at any given time.
 * As aligners typically use all available cores, throughput can be improved by
 * ensuring only one instance is running at any given time.
 */
public class SynchronisedFastqAligner implements FastqAligner {
	private final FastqAligner aligner;
	public SynchronisedFastqAligner(FastqAligner aligner) {
		this.aligner = aligner;
	}
	public synchronized void align(File fastq, File output, File reference, int threads) throws IOException {
		aligner.align(fastq, output, reference, threads);
	}
}
