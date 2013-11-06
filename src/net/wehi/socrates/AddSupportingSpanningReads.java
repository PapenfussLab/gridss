import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.StringTokenizer;

import javax.swing.text.html.HTMLDocument.HTMLReader.IsindexAction;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileReader.ValidationStringency;



public class AddSupportingSpanningReads {
	
	private static boolean is5primeSoftClip(SAMRecord s){
		int sc_start = s.getAlignmentStart() - s.getUnclippedStart();
		if (sc_start > 0)
			return true;
		return false;
	}
	private static int softclipLength5prime(SAMRecord s){
		int sc_start = s.getAlignmentStart() - s.getUnclippedStart();
		return sc_start;
	}
	private static int softclipLength3prime(SAMRecord s){
		int sc_end = s.getUnclippedEnd() - s.getAlignmentEnd();
		return sc_end;
	}
	private static boolean is3primeSoftClip(SAMRecord s){
		int sc_end = s.getUnclippedEnd() - s.getAlignmentEnd();
		if (sc_end > 0)
			return true;
		return false;
	}
	private static boolean isSoftClip(SAMRecord s){
		if (s.getUnclippedEnd() != s.getAlignmentEnd() || s.getUnclippedStart() != s.getAlignmentStart())
			return true;
		return false;
	}
	private static boolean isInteresingSoftclip(SAMRecord s, int breakpoint){
		if(s.getAlignmentStart()-2 <= breakpoint && s.getAlignmentEnd()+2>=breakpoint && (softclipLength3prime(s) > 4 || softclipLength5prime(s) > 4))
			return true;
//		if(s.getAlignmentStart() breakpoint  && )
//			return true;
		return false;
	}
	
	private static boolean isAnomalousPair(SAMRecord s, int lowerBound, int upperBound, int breakpoint){
		
		if (!s.getProperPairFlag())
			return true;
		
		if(!s.getReferenceName().equals(s.getMateReferenceName()) )
			return true;
		
		int iSize = Math.abs(s.getInferredInsertSize());
		if(iSize < lowerBound || iSize > upperBound){
			if(s.getAlignmentStart() < breakpoint && s.getInferredInsertSize() > 0)
				return true;
			if(s.getAlignmentEnd() > breakpoint && s.getInferredInsertSize() < 0)
				return true;
		}
		
		return false;
	}
	
        private static boolean mateMatch(SAMRecord s, String chr, int start, String dir, int upper, int readToBreakpointDist){
		
		
                int end;

		if(dir.equals("-"))
		{
                        end=start+1+(upper-readToBreakpointDist);
		        start=start+1;
                }
                else
                {			
		        end=start-1;
                        start=start-1-(upper-readToBreakpointDist);
                }
		if(s.getMateReferenceName().equals(chr) &&
                        s.getMateAlignmentStart()>=start &&
                        (s.getMateAlignmentStart()+100)<=end)
                {
                        return true;
                }
		return false;
	}
	
	private static boolean isConcordantPair(SAMRecord s, int lowerBound, int upperBound, int breakpoint){
		if (!s.getProperPairFlag())
			return false;
		if(!s.getReferenceName().equals(s.getMateReferenceName()) )
			return false;
		int iSize = Math.abs(s.getInferredInsertSize());
		if(iSize >= lowerBound && iSize <= upperBound){
			if(s.getAlignmentStart() < breakpoint && s.getInferredInsertSize() > 0)
				return true;
			if(s.getAlignmentEnd() > breakpoint && s.getInferredInsertSize() < 0)
				return true;
		}
		return false;
	}
	
	private static boolean isAlignedAcrossBreakpoint(SAMRecord s, int breakpointPosition){
		if(s.getAlignmentStart() < breakpointPosition-4 && s.getAlignmentEnd() > breakpointPosition+4)
			return true;
		return false;
	}
	

	private static int fetchStuff(
			String bp1chr, int bp1start, String bp1dir, String bp2chr, int bp2start, String bp2dir,SAMFileReader bamFile, int lowerInsertBound, int upperInsertBound) {
		String chr;
                int start;
                int end;

		if(bp1dir.equals("-"))
		{
		        start=bp1start+1;
                        end=bp1start+1+(upperInsertBound-100);
                        chr=bp1chr;
                }
                else
                {			
		        end=bp1start-1;
                        start=bp1start-1-(upperInsertBound-100);
                        chr=bp1chr;
                }

		SAMRecordIterator iter = bamFile.queryContained(chr, start, end);
		
		
		int spanningReadCount=0;
		
                int readToBreakpointDist;
		for(SAMRecord s; iter.hasNext();){
			s = iter.next();
                        if(bp1dir.equals("-"))
                        {
                                readToBreakpointDist=s.getAlignmentStart()-bp1start;
                        }
                        else
                        {
                                readToBreakpointDist=bp1start-s.getAlignmentEnd();
                        }

                        if(mateMatch(s,bp2chr,bp2start,bp2dir,upperInsertBound,readToBreakpointDist))
                        {
                                spanningReadCount++;
                        }
		}
		iter.close();
		return spanningReadCount;
	}
	
	private static String resultsToString(int[] result){
		if (result.length==5)
			return result[0]+"\t"+result[1]+"\t"+result[2]+"\t"+result[3]+"\t"+result[4];
		return result[0]+"\t"+result[1]+"\t"+result[2];
	}
	
	public static void main(String[] args) throws IOException {
		
		if (args.length < 3) {
			System.err.println("Usage: <bam file> <mean fragment length> <SD fragment length> <Socrates output>");
			System.exit(0);
		}
		
		double mean = Double.parseDouble(args[1]);
		double std = Double.parseDouble(args[2]);
		int upper = (int) Math.round(mean + 3*std);
		int lower = (int) Math.round(mean - 3*std);
		
		
		String bamfilename = args[0];
		File bamfile = new File(bamfilename);
		//System.out.println(bamfilename);	
		final SAMFileReader reader = new SAMFileReader(bamfile);
	        reader.setValidationStringency(ValidationStringency.SILENT);	
		
		BufferedReader vars = new BufferedReader(new FileReader(args[3]));
		String line;
		
		while( (line = vars.readLine()) != null) {
                        //System.out.println("*"+line);
			String[] cols = line.split("\t");
			StringTokenizer bp1 = new StringTokenizer(cols[3],":");
			String bp1dir = cols[4];
			StringTokenizer bp2 = new StringTokenizer(cols[15],":");
			String bp2dir = cols[16];

			String bp1chr = bp1.nextToken();
			int bp1start = Integer.parseInt(bp1.nextToken());
			int bp1end = bp1start;
			String bp2chr = bp2.nextToken();
			int bp2start = Integer.parseInt(bp2.nextToken());
			int bp2end = bp1start;
			
			//System.out.println(bp1chr+":"+bp1start+"\t"+bp1dir+"\t"+bp2chr+":"+bp2start+"\t"+bp2dir);
			
			int results = fetchStuff(bp1chr, bp1start, bp1dir, bp2chr, bp2start, bp2dir, reader, lower, upper);
			System.out.println(line+"\t"+results);
		}
		
	}
	
}
