package htsjdk.samtools.fastq;

import java.io.File;
import java.io.Flushable;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.util.IOUtil;

/**
 * BasicFastqWriter that does not flush on very write() call.
 *
 */
public class NonFlushingBasicFastqWriter implements FastqWriter,Flushable {
	// TODO refactor into FastqConstants
	private static final Charset FASTQ_ENCODING = StandardCharsets.UTF_8;
	// TODO refactor into FastqConstants
	public static final byte[] SEQUENCE_HEADER = FastqConstants.SEQUENCE_HEADER.getBytes(FASTQ_ENCODING);
    public static final byte[] QUALITY_HEADER =  FastqConstants.QUALITY_HEADER.getBytes(FASTQ_ENCODING);
    public static final byte[] FIRST_OF_PAIR =  FastqConstants.FIRST_OF_PAIR.getBytes(FASTQ_ENCODING);
    public static final byte[] SECOND_OF_PAIR =  FastqConstants.SECOND_OF_PAIR.getBytes(FASTQ_ENCODING);
    // Note: this differs from the htsjdk implementation which writes EOL using the system default
    // This seems preferable to having some tools break if we write using CRLF 
    private static final byte[] FASTQ_EOL = new byte[] { '\n' };
    //private static final byte[] FASTQ_EOL = System.lineSeparator().getBytes(FASTQ_ENCODING);
    
	private final String path;
	private final OutputStream stream;

	public NonFlushingBasicFastqWriter(final File file) {
		this(file, IOUtil.maybeBufferOutputStream(IOUtil.openFileForWriting(file)));
	}

	private NonFlushingBasicFastqWriter(final File file, final OutputStream writer) {
		this.path = (file != null? file.getAbsolutePath(): "");
		this.stream = writer;
	}

	public NonFlushingBasicFastqWriter(final OutputStream writer) {
		this(null, writer);
	}

	@Override
	public void write(final FastqRecord rec) {
		try {
			writeStream(stream, rec);
		} catch (IOException e) {
			throw new SAMException("Error in writing fastq file " + path, e);
		}
	}
	// TODO: refactor into FastqEncoder.write()
	private static void writeStream(final OutputStream stream, final FastqRecord record) throws IOException {
		final String readName = record.getReadName();
		final String readString = record.getReadString();
		final String qualHeader = record.getBaseQualityHeader();
		final String qualityString = record.getBaseQualityString();
		stream.write(SEQUENCE_HEADER);
		if (readName != null) {
			stream.write(readName.getBytes(FASTQ_ENCODING));
		}
		stream.write(FASTQ_EOL);
		if (readString != null) {
			stream.write(readString.getBytes(FASTQ_ENCODING));
		}
		stream.write(FASTQ_EOL);
		stream.write(QUALITY_HEADER);
		if (qualHeader != null) {
			stream.write(qualHeader.getBytes(FASTQ_ENCODING));
		}
		stream.write(FASTQ_EOL);
		if (qualityString != null) {
			stream.write(qualityString.getBytes(FASTQ_ENCODING));
		}
		stream.write(FASTQ_EOL);
	}
    @Override
    public void flush() throws IOException {
    	// Design decision: should we wrap this with a SAMException like we do in write() and close()?
    	stream.flush();
    }

    @Override
    public void close() {
    	try {
    		flush();
        	stream.close();
		} catch (IOException e) {
			throw new SAMException("Error in writing fastq file " + path, e);
		}
    }
}
