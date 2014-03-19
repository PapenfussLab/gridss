/**
 * 
 */
package au.edu.wehi.socrates.util;

//import net.sf.picard.*;

import net.sf.samtools.*;
import net.wehi.socrates.SOCRATES;

import java.awt.geom.Point2D;
import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.ArrayList;


/**
 * @author hsu
 * @author schroeder
 * Created on Jan 18, 2013
 */
@Deprecated() // DC: use net.sf.picard.metrics
public class SAMFileInfo {
	public ArrayList<String> libraries = new ArrayList<String>();
	public ArrayList<String> sequenceNames = new ArrayList<String>();
	private HashMap<String,Integer> sequenceName2ID = new HashMap<String,Integer>();
	private HashMap<Integer,String> sequenceID2Name = new HashMap<Integer,String>();
	public HashMap<String,Integer> sequenceLengths = new HashMap<String,Integer>();
	
	public SAMFileInfo(SAMFileReader reader) {
		SAMFileHeader header = reader.getFileHeader();
		
		// add library IDs
		for (SAMReadGroupRecord rg : header.getReadGroups()) {
			libraries.add( rg.getReadGroupId() );
		}
		if (libraries.size()==0) libraries.add( "null" );
		
		// add sequence name to id
		for (SAMSequenceRecord s : header.getSequenceDictionary().getSequences()) {
			String sname = s.getSequenceName();
			sequenceNames.add(sname);
			sequenceLengths.put( sname, s.getSequenceLength() );
			sequenceName2ID.put(sname, s.getSequenceIndex());
			sequenceID2Name.put(s.getSequenceIndex(), sname);
		}
	}
	
	public SAMFileInfo(String readerFilename) {
		SAMFileReader reader = new SAMFileReader(new File(readerFilename));
		SAMFileHeader header = reader.getFileHeader();
		
		// add library IDs
		for (SAMReadGroupRecord rg : header.getReadGroups()) {
			libraries.add( rg.getReadGroupId() );
		}
		if (libraries.size()==0) libraries.add( "null" );
		
		// add sequence name to id
		for (SAMSequenceRecord s : header.getSequenceDictionary().getSequences()) {
			String sname = s.getSequenceName();
			sequenceNames.add(sname);
			sequenceLengths.put( sname, s.getSequenceLength() );
			sequenceName2ID.put(sname, s.getSequenceIndex());
			sequenceID2Name.put(s.getSequenceIndex(), sname);
		}
		reader.close();
	}
	
	public int getSequenceLength(int seqIndex) {
		return sequenceLengths.get(getSequenceName(seqIndex));
	}
	
	public int getSequenceLength(String seqName) {
		return sequenceLengths.get(seqName);
	}

	public int getSequenceIndex(String seqName) {
		if (sequenceName2ID.containsKey(seqName)) return sequenceName2ID.get(seqName);
		return -1;
	}
	
	public String getSequenceName(int id) {
		if (sequenceID2Name.containsKey(id)) return sequenceID2Name.get(id);
		return null;
	}

	public static HashMap<String, Integer> getLibraryReadLengths(String metricFile) {
		try {
			HashMap<String, Integer> readLengths = new HashMap<String, Integer>();
			java.io.BufferedReader reader = new java.io.BufferedReader(new java.io.FileReader(metricFile));
			
			String line;
			while ( (line=reader.readLine()) != null ) {
				String[] tokens = line.split("\t", -1);
				String rg = tokens[0];
				int l = Integer.parseInt(tokens[tokens.length-1]);
				if (SOCRATES.verbose) System.out.println("Read length for read group " + rg + " = " + l);
			}
			
			reader.close();
			return readLengths;
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}
	
	public static HashMap<String, Point2D.Float> getInsertMeanStdev(String metricFile) {
		try {
			HashMap<String, Point2D.Float> isizeBounds = new HashMap<String, Point2D.Float>();
			java.io.BufferedReader reader = new java.io.BufferedReader(new java.io.FileReader(metricFile));
			
			String line;
			while ( (line=reader.readLine()) != null ) {
				String[] tokens = line.split("\t", -1);
				String rg = tokens[0];
				
				float mean = tokens[1].equals("null") ? Float.NaN : Float.parseFloat(tokens[1]);
				float stdev = tokens[2].equals("null") ? Float.NaN : Float.parseFloat(tokens[2]);
				isizeBounds.put( rg, new Point2D.Float(mean,stdev) );
				if (SOCRATES.verbose) System.out.println("Insert size bounds for read group " + rg + ": mean = " + mean + "\tstdev = " + stdev);
			}
			
			reader.close();
			return isizeBounds;
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}

	public static float[] getInsertMeanStdReadlen(String metricsFile) {
		try{
			java.io.BufferedReader reader = new java.io.BufferedReader(new java.io.FileReader(metricsFile));
			String line;
			float sumIsize = 0;
			float sumStd = 0;
			float sumReadlen = 0;
			int count = 0;
                        while ( (line=reader.readLine()) != null ) {
				String[] tokens = line.split("\t", -1);
                                String rg = tokens[0];

                                float mean = tokens[1].equals("null") ? Float.NaN : Float.parseFloat(tokens[1]);
                                float stdev = tokens[2].equals("null") ? Float.NaN : Float.parseFloat(tokens[2]);
				float readLen = tokens[3].equals("null") ? Float.NaN : Float.parseFloat(tokens[3]);
	
				sumIsize += mean;
				sumStd += stdev;
				sumReadlen += readLen;
				count ++;
	
				if(! (mean > sumIsize/count - sumStd/count && mean < sumIsize/count + sumStd/count) ) {
					System.err.println("Read groups do not have homogeneous insert sizes. Filtering could be unreliable with mixed data.	");			
				}
				if(  readLen != sumReadlen/count){
					System.err.println("Read groups have different read lengths. Taking average. Filtering could be unreliable with mixed data");
				}
			}
			reader.close();
			return new float[] {sumIsize/count, sumStd/count, sumReadlen/count};
		} catch (Exception e) {
                        e.printStackTrace();
                        return null;
                }
        }

					
	public static void writeBAMMetrics(SAMFileReader reader, ArrayList<String> libraries, String outputFile, int sampleSize) {
		if (SOCRATES.verbose) System.out.println("Writing metrics file");
		HashMap<String, Point2D.Float> isizeMeanStdev = getInsertSizeMeanStdDev(reader, libraries, sampleSize);
		HashMap<String, Integer> readLengths = getLibraryReadLengths(reader, libraries);
		
		try {
			File ofile = new File(outputFile);
			java.io.FileWriter out = new java.io.FileWriter( ofile.getName() );
			java.io.PrintWriter p = new java.io.PrintWriter( out );
			
			for (String rgId : isizeMeanStdev.keySet()) {
				Point2D.Float bounds = isizeMeanStdev.get(rgId);
				Integer length = readLengths.get(rgId);
				if (bounds!=null) p.println(rgId + "\t" + bounds.x + "\t" + bounds.y + "\t" + length);
				else p.println(rgId + "\tnull\tnull\t" + length);
			}
			p.close();
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private static HashMap<String, Point2D.Float> getInsertSizeMeanStdDev(SAMFileReader reader, ArrayList<String> libraries, int sampleSize) {
		HashMap<String,ArrayList<Integer>> isizes = new HashMap<String,ArrayList<Integer>>();
		HashSet<String> completed = new HashSet<String>();

		// store absolute insert sizes for each library
		SAMRecordIterator iter = reader.iterator();
		while (iter.hasNext()) {
			SAMRecord aln = iter.next();
			SAMReadGroupRecord rgr = aln.getReadGroup();
			String rg = (rgr==null) ? "null" : rgr.getReadGroupId();

			if (aln.getReadPairedFlag()) {
				if (!aln.getProperPairFlag()) continue;
				ArrayList<Integer> is = isizes.get(rg);
				
				if (is==null) {
					isizes.put( rg, new ArrayList<Integer>() );
					is = isizes.get(rg);
				}
				
				if (is.size() >= sampleSize) continue;
				
				is.add( Math.abs(aln.getInferredInsertSize()) );
				
				if (is.size() == sampleSize) {
					completed.add( rg );
					if (completed.size()==libraries.size()) break;
				}
			} else {
				completed.add( rg );
				if (completed.size()==libraries.size()) break;
			}
		}
		iter.close();
		
		// compute mean insert size
		HashMap<String, Point2D.Float> isizeMeanStd = new HashMap<String, Point2D.Float>();
		
		long[] sum = new long[libraries.size()];
		for (int i=0; i<libraries.size(); i++) {
			String rg = libraries.get(i);
			if (isizes.get(rg)!=null) 
				for (Integer v : isizes.get(rg)) sum[i] += v;
		}

		float[] sum_dev = new float[libraries.size()];
		// compute standard deviation of insert size
		for (int i=0; i<libraries.size(); i++) {
			String rg = libraries.get(i);
			ArrayList<Integer> is = isizes.get(rg);
			
			if (is!=null) {
				float avg = sum[i]/(float)(is.size());
				for (Integer v : is) {
					float dev = v - avg;
					sum_dev[i] += dev*dev;
				}
				float stdev = (float)(Math.sqrt(sum_dev[i]/((float)is.size()-1) ));
				isizeMeanStd.put(rg, new Point2D.Float(avg, stdev));
				if (SOCRATES.verbose) System.out.println("Insert size bounds for read group " + rg + ": mean = " + avg + "\tstdev = " + stdev);
			} else {
				isizeMeanStd.put(rg, null);
				if (SOCRATES.verbose) System.out.println("Insert size bounds for read group " + rg + ": mean = " + "null" + "\tstdev = " + "null");
			}
		}

		return isizeMeanStd;
	}
	
	private static HashMap<String, Integer> getLibraryReadLengths(SAMFileReader reader, ArrayList<String> libraries) {
		HashMap<String, Integer> readLengths = new HashMap<String, Integer>();
		
		SAMRecordIterator iter = reader.iterator();
		while (iter.hasNext()) {
			SAMRecord aln = iter.next();
			
			SAMReadGroupRecord rgr = aln.getReadGroup();
			String rg = rgr==null ? "null" : rgr.getReadGroupId();
			if (readLengths.containsKey(rg)) continue;
			readLengths.put( rg, aln.getReadLength() );
			if (SOCRATES.verbose) System.out.println("Read length for read group " + rg + " = " + aln.getReadLength());
			if (readLengths.size() == libraries.size()) break;
		}
		iter.close();
		
		return readLengths;
	}
	
	public static void main(String[] args) {
		if (args.length<1) return;
		//final File bamFile = new File(args[0]);
		//final SAMFileReader reader = new SAMFileReader(bamFile);
		//SAMFileInfo info = new SAMFileInfo(reader);
		//writeBAMMetrics(reader, info.libraries, bamFile.getName()+".metrics", 10000);
		float[] f = getInsertMeanStdReadlen(args[0]);
		System.out.println(f[0]+"\t"+f[1]+"\t"+f[2]);
	}
}
