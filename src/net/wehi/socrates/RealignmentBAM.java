/**
 * 
 */
package net.wehi.socrates;

import java.io.File;
//import java.io.FileInputStream;
//import java.io.BufferedReader;
//import java.io.InputStreamReader;
//import java.util.zip.GZIPInputStream;
//import java.util.HashMap;


import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.BAMRecordCodec;
import net.sf.samtools.SAMRecordCoordinateComparator;
import net.sf.samtools.util.SortingCollection;
import net.sf.samtools.util.CloseableIterator;
import net.wehi.socrates.util.ProgressVerbose;
//import net.wehi.socrates.util.Utilities;

/**
 * @author hsu
 *
 * Created on Jan 23, 2013
 */
public class RealignmentBAM {
	public RealignmentBAM() {}

	
	public static void makeRealignmentBAM(String inputBAMFilename, String outputBamFilename) {
		final SAMFileReader sam = (inputBAMFilename.equals("-")) ? new SAMFileReader( System.in ) : new SAMFileReader( new File(inputBAMFilename) );
		SAMFileReader.setDefaultValidationStringency(SAMFileReader.ValidationStringency.SILENT);
		
		ProgressVerbose p = new ProgressVerbose("Sorting alignments", 1000000, SOCRATES.verbose);
		try {
			SortingCollection<SAMRecord> buffer = SortingCollection.newInstance(SAMRecord.class, new BAMRecordCodec(sam.getFileHeader()), 
																				new SAMRecordCoordinateComparator(), 500000, new File("."));
			
			SAMRecordIterator iter = sam.iterator();
			while (iter.hasNext()) {
				SAMRecord aln = iter.next();
				
				if (aln.getNotPrimaryAlignmentFlag()) continue;
				
				p.stepProgress(" alignments added to sorting collection");
				String readName = aln.getReadName();
				String[] tokens = readName.split("&");
				assert tokens.length == 5;
				
				aln.setReadName( tokens[0] );
				aln.setFirstOfPairFlag(true);
				aln.setProperPairFlag(false);
				aln.setReadPairedFlag(true);
				aln.setMateReferenceIndex( Integer.parseInt(tokens[1]) );
				aln.setMateAlignmentStart( Integer.parseInt(tokens[2]) );
				aln.setMateNegativeStrandFlag( tokens[3].equals("-") );
				aln.setMateUnmappedFlag( !(tokens[4].equals("1")) );
				
				aln.setAttribute("ZS", tokens[5]);
				buffer.add(aln);
			}
			iter.close();
			p.end(" alignments added to sorting collection");
			
			ProgressVerbose p2 = new ProgressVerbose("Creating realignment BAM file", 1000000, SOCRATES.verbose);
			SAMFileWriterFactory.setDefaultCreateIndexWhileWriting(true);
			sam.getFileHeader().setSortOrder(SortOrder.coordinate);
			SAMFileWriter bam = new SAMFileWriterFactory().makeBAMWriter( sam.getFileHeader(), true, new File(outputBamFilename) );
			
			CloseableIterator<SAMRecord> iter2 = buffer.iterator();
			while (iter2.hasNext()) {
				SAMRecord aln = iter2.next();
				bam.addAlignment(aln);
				p2.stepProgress(" alignments written");
			}
			p2.end(" alignments written");
			sam.close();
			bam.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}


	class Info {
		public int chromIdx, pos;
		public char orientation;
		public byte[] seq;
		public boolean ideal;
		
		public Info() {};
		
		public Info(int chromIdx, int pos, char orientation, byte[] seq, char ideal) {
			this.chromIdx = chromIdx;
			this.pos = pos;
			this.orientation = orientation;
			this.seq = seq;
			this.ideal = (ideal=='1');
		}
	}

    public static void printHelp() {
        System.err.println("usage: RealignmentBAM input_bam output_bam");
        System.err.println(" set input_bam to \"-\" to accept input from stdin");
    }

	
	public static void main(String[] args) {
		if (args.length==2) {
            SOCRATES.verbose = false;
            System.err.println("\nAdd anchor information into re-alignment BAM file");
            System.err.println("  input BAM file:\t" + args[0] );
            System.err.println("  output BAM file:\t" + args[1] );
			makeRealignmentBAM(args[0], args[1]);
		}
		else {
			for (int i=0; i<args.length; i++)
				System.err.println(args[i]);
            printHelp();
            System.exit(1);
		}
	}
}

