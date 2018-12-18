package au.edu.wehi.idsv;

import au.edu.wehi.idsv.bed.IntervalBed;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import au.edu.wehi.idsv.util.FileHelper;
import htsjdk.samtools.*;

import java.io.File;
import java.io.IOException;

public class LinearScanFullReadExtractor extends FullReadExtractor {
    public LinearScanFullReadExtractor(LinearGenomicCoordinate lgc, IntervalBed bed, boolean extractMates, boolean extractSplits) {
        super(lgc, bed, extractMates, extractSplits);
    }

    @Override
    public void extract(File input, File output, int workerThreads) throws IOException {
        File tmpOut = gridss.Defaults.OUTPUT_TO_TEMP_FILE ? FileSystemContext.getWorkingFileFor(output) : output;
        try (SamReader reader = SamReaderFactory.makeDefault().open(input)) {
            SAMFileHeader header = reader.getFileHeader();
            try (AsyncBufferedIterator<SAMRecord> asyncIt = new AsyncBufferedIterator<>(reader.iterator(), input.getName())) {
                try (SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, tmpOut)) {
                    while (asyncIt.hasNext()) {
                        SAMRecord r = asyncIt.next();
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
