/**
 * 
 */
package au.edu.wehi.socrates.util;

import java.io.File;
import java.io.IOException;

import net.sf.samtools.seekablestream.SeekableStream;



/**
 * @author hsu
 *
 * Created on Jan 25, 2013
 */
public class SeekableMemoryStream extends SeekableStream {
	protected MemoryMappedFile mem;
	protected long pos;
	
	public SeekableMemoryStream(MemoryMappedFile file) {
		mem = file;
		pos = 0;
	}

	/* (non-Javadoc)
	 * @see net.sf.samtools.util.SeekableStream#close()
	 */
	@Override
	public void close() throws IOException {
	}

	/* (non-Javadoc)
	 * @see net.sf.samtools.util.SeekableStream#eof()
	 */
	@Override
	public boolean eof() throws IOException {
		return (pos>mem.length());
	}

	/* (non-Javadoc)
	 * @see net.sf.samtools.util.SeekableStream#getSource()
	 */
	@Override
	public String getSource() {
		return "Memory mapped " + mem.getFilename();
	}
	
	public File getSourceFile() { return mem.getFile(); }

	/* (non-Javadoc)
	 * @see net.sf.samtools.util.SeekableStream#length()
	 */
	@Override
	public long length() {
		return mem.length();
	}

	/* (non-Javadoc)
	 * @see net.sf.samtools.util.SeekableStream#read(byte[], int, int)
	 */
	@Override
	public int read(byte[] buffer, int offset, int len) throws IOException {
		if (pos>mem.length()) return -1;
		
		int read = 0;
		int startPage = (int)(pos/MemoryMappedFile.PAGE_LIMIT);
		int endPage = (int)((pos+len)/MemoryMappedFile.PAGE_LIMIT);
		
		if (startPage==endPage) { // all data on the same page of memory
			int idx  = (int)(pos%MemoryMappedFile.PAGE_LIMIT);
			byte[] mappedPage = mem.getMappedPage(startPage);
			int max = offset+len;
			for (int p=offset; p<max; p++) {
//				try {
				buffer[p] = mappedPage[idx];
//				} catch (Exception e) {
//					System.err.println("pages:\t" + mem.getPageCount() );
//					System.err.println("page:\t" + startPage);
//					System.err.println("page length:\t" + mappedPage.length );
//					System.err.println("page index:\t" + idx);
//					System.err.println("buffer length:\t" + buffer.length );
//					System.err.println("buffer index:\t" + p);
//					throw new IOException("ERROR");
//				}
				idx++;
				pos++;
				read++;
				if (pos>mem.length()) break; 
			}
		} else { // data spans multiple pages of memory
			int oldPage = -1;
			byte[] mappedPage = null;
			for (int p=offset; p<offset+len; p++) {
				int page = (int)(pos/MemoryMappedFile.PAGE_LIMIT);
				if (page!=oldPage) {
					oldPage = page;
					mappedPage = mem.getMappedPage(page);
				}
				buffer[p] = mappedPage[(int)(pos%MemoryMappedFile.PAGE_LIMIT)];
				pos++;
				read++;
				if (pos>mem.length()) break; 
			}
		}
		return read;
	}

	/* (non-Javadoc)
	 * @see net.sf.samtools.util.SeekableStream#seek(long)
	 */
	@Override
	public void seek(long p) throws IOException {
		if (p>mem.length()) {
			throw new IOException("Exceed file size");
		}
		pos = p;
	}

	/* (non-Javadoc)
	 * @see java.io.InputStream#read()
	 */
	@Override
	public int read() throws IOException {
		if (pos>mem.length()) {
			return -1;
		}
		
		int ret = mem.getMappedPage((int)(pos/MemoryMappedFile.PAGE_LIMIT))[(int)(pos%MemoryMappedFile.PAGE_LIMIT)];
		pos++;
		
		return ret;
	}
	@Override
	public long position() throws IOException {
		return pos;
	}
}
