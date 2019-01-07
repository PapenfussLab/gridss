package au.edu.wehi.idsv;

import au.edu.wehi.idsv.bed.IntervalBed;
import htsjdk.samtools.*;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;

public class IndexedLookupFullReadExtractor extends FullReadExtractor {
    private final int regionPaddingSize;
    public IndexedLookupFullReadExtractor(LinearGenomicCoordinate lgc, IntervalBed bed, boolean extractMates, boolean extractSplits, int regionPaddingSize) {
        super(lgc, bed, extractMates, extractSplits);
        this.regionPaddingSize = regionPaddingSize;
    }
    public void extract(File input, File output, int workerThreads) throws IOException {
        IntervalBed lookupIntervals = new IntervalBed(getLinearGenomicCoordinate().getDictionary(), getLinearGenomicCoordinate());
        for (QueryInterval qi : getRegionBed().asQueryInterval()) {
            lookupIntervals.addInterval(expandBy(qi, regionPaddingSize, getLinearGenomicCoordinate().getDictionary()));
        }
        // todo: convert QueryInterval to bam chunk
        // per worker parallel
        // proper bam merge, or concat?
        // while list not empty
            // process whole list
                // add mate/split locations to new list
            // list = new list
        BAMFileReader r;
        throw new RuntimeException("NYI");
    }
    private IntervalBed processInterval(QueryInterval qi, File input, BlockingQueue<Collection<SAMRecord>> outputBuffer) {
        throw new RuntimeException("NYI");
    }
    public List<Chunk> getChunks(File input) {
        throw new RuntimeException("NYI");

    }
    private static QueryInterval expandBy(QueryInterval qi, int regionPaddingSize, SAMSequenceDictionary dict) {
        return new QueryInterval(
                qi.referenceIndex,
                Math.max(0, qi.start),
                Math.min(dict.getSequence(qi.referenceIndex).getSequenceLength(), qi.end));
    }
}
