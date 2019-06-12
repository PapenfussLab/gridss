package au.edu.wehi.idsv;

import au.edu.wehi.idsv.bed.IntervalBed;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import au.edu.wehi.idsv.util.FileHelper;
import htsjdk.samtools.*;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;

import java.io.File;
import java.io.IOException;

public class FullReadExtractor extends ReadExtractor {
    private static final Log log = Log.getInstance(FullReadExtractor.class);

    public FullReadExtractor(LinearGenomicCoordinate lgc, IntervalBed bed, boolean extractMates, boolean extractSplits) {
        super(lgc, bed, extractMates, extractSplits);
    }

    @Override
    public void extract(File input, File output, int workerThreads) throws IOException {
        File tmpOut = gridss.Defaults.OUTPUT_TO_TEMP_FILE ? FileSystemContext.getWorkingFileFor(output) : output;
        try (SamReader reader = SamReaderFactory.makeDefault().open(input)) {
            SAMFileHeader header = reader.getFileHeader();
            try (AsyncBufferedIterator<SAMRecord> asyncIt = new AsyncBufferedIterator<>(reader.iterator(), input.getName())) {
                ProgressLoggingSAMRecordIterator it = new ProgressLoggingSAMRecordIterator(asyncIt, new ProgressLogger(log));
                try (SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, tmpOut)) {
                    while (it.hasNext()) {
                        SAMRecord r = it.next();
                        if (shouldExtract(r)) {
                            writer.addAlignment(r);
                        }
                    }
                }
            }
        }
        if (tmpOut != output) {
            FileHelper.move(tmpOut, output, true);
        }
    }
}
