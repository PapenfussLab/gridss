package au.edu.wehi.idsv.alignment;

import au.edu.wehi.idsv.picard.ReferenceLookup;
import htsjdk.samtools.SAMSequenceDictionary;

import java.io.File;
import java.io.IOException;

public interface FastqAligner {
	void align(File fastq, File output, File reference, int threads, SAMSequenceDictionary dict) throws IOException;
}
