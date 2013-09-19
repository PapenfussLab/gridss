/**
 * 
 */
package net.wehi.socrates.util;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.concurrent.atomic.AtomicReferenceArray;

import net.wehi.socrates.SOCRATES;

/**
 * @author hsu
 *
 * Created on Jan 31, 2013
 */
public class MemoryMappedFile {
	public static final int PAGE_LIMIT = 50*1024*1024; // using 50MB blocks
	protected File source;
	protected AtomicReferenceArray<byte[]> data;
	protected byte[][] eagerData;
	protected long fileSize;
	
	/**
	 * 
	 */
	public MemoryMappedFile(File file, boolean eagerLoad) {
		if (SOCRATES.verbose) {
			System.out.println("Mapping file to memory " + file.getName());
		}
		try {
			source = file;
			RandomAccessFile randomAccess = new RandomAccessFile(file, "r");
			fileSize = randomAccess.length();
			
			int pages = (int)(fileSize/PAGE_LIMIT)+1;
			eagerData = eagerLoad ? new byte[pages][] : null;
			
			if (eagerLoad) {
				for (int page=0; page<pages; page++) {
					int pageSize = (page!=pages-1) ? PAGE_LIMIT : (int)(fileSize%PAGE_LIMIT);
					eagerData[page] = new byte[pageSize];
					randomAccess.readFully(eagerData[page]);
				}
				data = null;
			} else {
				data = new AtomicReferenceArray<byte[]>(pages);
			}
			randomAccess.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		if (SOCRATES.verbose) {
			System.out.println("Mapping file to memory - DONE");
		}
	}
	
//	private void readAll() throws IOException {
//		RandomAccessFile randomAccess = new RandomAccessFile(source, "r");
//		for (int page=0; page<data.length(); page++) {
//			if (page!=data.length()-1) readData(page*PAGE_LIMIT, PAGE_LIMIT);
//			else readData(page*PAGE_LIMIT, (int)(fileSize%PAGE_LIMIT));
//		}
//		randomAccess.close();
//	}

	private byte[] readData(long pos, int size) {
		byte[] b = new byte[size];
		try {
			RandomAccessFile randomAccess = new RandomAccessFile(source, "r");
			randomAccess.seek(pos);
			randomAccess.readFully(b);
			randomAccess.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return b;
	}
	

	public byte[] getMappedPage(int page) throws IOException {
		if (eagerData != null) return eagerData[page];
		
		data.compareAndSet(page, null, (page!=data.length()-1) ? readData(page*PAGE_LIMIT, PAGE_LIMIT) :
																 readData(page*PAGE_LIMIT, (int)(fileSize%PAGE_LIMIT)) );
		return data.get(page);
	}
	
	public String getFilename() { return source.getName(); }
	
	public File getFile() { return source; }
	
	public long length() { return fileSize; }
}
