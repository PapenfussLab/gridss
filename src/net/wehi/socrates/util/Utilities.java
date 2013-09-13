/**
 * 
 */
package net.wehi.socrates.util;

import java.io.File;
import java.util.Arrays;

/**
 * @author hsu
 *
 * Created on Jan 18, 2013
 */

public class Utilities {
	/**
	 * Ensures the given file exists in the file system.
	 * @param filename
	 * @return
	 */
	public static File assertExists(String filename) {
		File f = new File(filename);
		if (!f.exists()) {
			System.err.println("Required file "+f.getName()+" does not exist. Exiting...");
			System.exit(1);
		}
		return f;
	}
	
	/**
	 * Deletes BAM file and its index file if exists.
	 * @param bamFile
	 */
	public static void deleteBAM(File bamFile) {
		String fn = bamFile.getAbsolutePath();
		File indexFile = new File( fn.substring(0,fn.lastIndexOf('.'))+".bai" );
		if (indexFile.exists()) indexFile.delete();
		indexFile = new File( fn+".bai" );
		if (indexFile.exists()) indexFile.delete();
		if (bamFile.exists()) bamFile.delete();
	}
	
	
	public static char nucleotideComplement(char n) {
		switch (n) {
		case 'A': return 'T';
		case 'T': return 'A';
		case 'C': return 'G';
		case 'G': return 'C';
		case 'a': return 't';
		case 't': return 'a';
		case 'c': return 'g';
		case 'g': return 'c';
		default: return 'N';
		}
	}	

	public static int nucleotideToInt(char nt) {
		switch (nt) {
		case 'A': return 1;
		case 'T': return 2;
		case 'C': return 3;
		case 'G': return 4;
		default: return 0;
		}
	}
	
	public static char indexToNucleotide(int n) {
		switch (n) {
		case 1: return 'A';
		case 2: return 'T';
		case 3: return 'C';
		case 4: return 'G';
		default: return 'N';
		}
	}
	
	public static int byteToIndex(byte nt) {
		switch (nt) {
		case 65: return 1;
		case 84: return 2;
		case 67: return 3;
		case 71: return 4;
		default: return 0;
		}
	}
	
	public static int byteToIndexComplement(byte nt) {
		switch (nt) {
		case 65: return 2;
		case 84: return 1;
		case 67: return 4;
		case 71: return 3;
		default: return 0;
		}
	}
	
	public static byte byteToByteComplement(byte nt) {
		switch (nt) {
		case 65: return 84;
		case 84: return 65;
		case 67: return 71;
		case 71: return 67;
		default: return nt;
		}
	}
	
	public static byte[] getReversedArray(byte[] seq) {
		byte[] s = new byte[seq.length];
		int p = seq.length-1;
		for (int i=0; i<seq.length; i++) {
			s[p] = seq[i];
			p--;
		}
		return s;
	}
	
	public static byte[] getReversedComplementArray(byte[] seq) {
		byte[] s = new byte[seq.length];
		int p = seq.length-1;
		for (int i=0; i<seq.length; i++) {
			s[p] = byteToByteComplement(seq[i]);
			p--;
		}
		return s;
	}
	
	public static void reverse(byte[] seq) {
		for (int i=0; i < seq.length/2; i++) {
			byte tmp = seq[i];
			seq[i] = seq[seq.length-i-1];
			seq[seq.length-i-1] = tmp;
		}
	}
	
	public static void reverseComplement(byte[] seq) {
		for (int i=0; i < seq.length/2; i++) {
			byte tmp = byteToByteComplement(seq[i]);
			seq[i] = byteToByteComplement(seq[seq.length-i-1]);
			seq[seq.length-i-1] = tmp;
		}
		if (seq.length%2==1) seq[seq.length/2] = byteToByteComplement(seq[seq.length/2]);
	}
	
	public static String sequenceByteToString(byte[] seq, boolean reverse) {
		if (seq==null || seq.length==0) return new String();
		
		StringBuilder s = new StringBuilder();
		for (byte b : seq) {
			if (reverse) {
				s.append( Utilities.nucleotideComplement( (char)b ) );
			} else {
				s.append( (char)b );
			}
		}
		if (reverse) {
			s = s.reverse();
		}
		return s.toString();
	}
	
	public static String qualityByteToString(byte[] qual, boolean reverse) {
		if (qual==null) return null;
		
		StringBuilder s = new StringBuilder();
		for (byte b : qual) {
			s.append( (char)(b+33) );
		}
		if (reverse) s = s.reverse();
		return s.toString();
	}
	
	public static byte[] stringToByte(String seq) {
		try {
			return seq.getBytes("US-ASCII");
		} catch (Exception e) {
			byte[] b = new byte[seq.length()];
			for (int i=0; i<seq.length(); i++) b[i] = (byte)seq.charAt(i);
			return b;
		}
	}
	
	public static boolean consensusStartsWith(byte[] consensus, byte[] seq, float identity) {
		if (consensus==null || seq==null || consensus.length==0 || seq.length==0) return false;
		int allowedMismatches = (int)(seq.length * (1-identity));
		
		int p = 0;
		for (byte b : seq) {
			if (p>=consensus.length) break;
			if (consensus[p]!=78 /*N*/ && b != consensus[p]) {
				allowedMismatches--;
				if (allowedMismatches < 0) return false;
			}
			p++;
		}
		
		return true;
	}
	
	public static boolean consensusStartsWith(byte[] consensus, byte[] seq, float identity, int consensusSkip, int seqSkip) {
		if (consensus==null || seq==null || consensus.length==0 || seq.length==0) return false;
		int allowedMismatches = (int)((seq.length-seqSkip) * (1-identity));
		
		int p = 0;
		for (int i=seqSkip; i<seq.length; i++) {
			byte b = seq[i];
			if (p+consensusSkip>=consensus.length) break;
			if (consensus[p+consensusSkip]!=78 /*N*/ && b != consensus[p+consensusSkip]) {
				allowedMismatches--;
				if (allowedMismatches < 0) return false;
			}
			p++;
		}
		
		return true;
	}
	
	public static boolean consensusEqual(byte[] consensus1, byte[] consensus2, float identity) {
		if (consensus1==null || consensus2==null || consensus1.length==0 || consensus2.length==0 || consensus1.length!=consensus2.length) return false;
		int allowedMismatches = (int)(consensus2.length * (1-identity));
		
		int p = 0;
		for (byte b : consensus1) {
			if (b==78 /*N*/ || consensus2[p]==78 /*N*/) continue;
			if (b != consensus2[p]) {
				allowedMismatches--;
				if (allowedMismatches < 0) return false;
			}
			p++;
		}
		
		return true;
	}
	
	public static byte[] concatenateByteArrays(byte[] a1, byte[] a2) {
		byte[] a3 = Arrays.copyOf(a1, a1.length+a2.length);
		System.arraycopy(a2, 0, a3, a1.length, a2.length);
		return a3;
	}
}
