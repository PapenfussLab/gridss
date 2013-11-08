package net.wehi.socrates;
/**
 * 
 */

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintWriter;
import java.io.FileOutputStream;
//import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeSet;
import java.util.HashSet;

import net.wehi.socrates.util.SAMFileInfo;
//import net.wehi.socrates.util.SAMRecordSummary;
import net.wehi.socrates.util.TabixReader;
import net.wehi.socrates.util.Utilities;

import net.wehi.socrates.util.MemoryMappedFile;
import net.wehi.socrates.util.SeekableMemoryStream;

//import net.sf.samtools.SAMFileReader;
//import net.sf.samtools.SAMRecord;
//import net.sf.samtools.SAMRecordIterator;

import java.io.File;

import java.util.concurrent.*;
//import java.util.concurrent.atomic.AtomicInteger;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;


/**
 * @author hsu
 *
 * Created on Mar 14, 2013
 */
@SuppressWarnings("static-access")
public class AnnotatePairedClusters implements Callable<Integer> {
	private static final Options options = new Options();
	static {
		Option help = new Option( "h", "help", false, "print this message" );
		Option verbose = new Option( "v", "verbose", false, "be verbose of progress" );
/*		Option threads = OptionBuilder.withArgName( "threads" )
										.hasArg()
										.withDescription("Number of threads to use [default: 3]")
										.withType( Number.class )
										.withLongOpt( "threads" )
										.create( 't' );*/
		Option normal = OptionBuilder.withArgName("normal")
				.hasArg()
				.withDescription("Socrates paired breakpoint calls for normal sample")
				.withType( String.class )
				.withLongOpt("normal")
				.create('n');
		Option flank = OptionBuilder.withArgName("flank")
				.hasArg()
				.withDescription("Normal breakpoint within FLANK nt of tumour breakpoint is considered as the same [default=10(nt)]")
				.withType( Integer.class )
				.withLongOpt( "flank" )
				.create( 'f' );
		Option rmsk = OptionBuilder.withArgName("repeatmask")
				.hasArg()
				.withDescription("UCSC repeat masker track file")
				.withType( String.class )
				.withLongOpt( "repeatmask" )
				.create( 'r' );
		
		//options.addOption( threads );
		options.addOption(normal);
		options.addOption(rmsk);
		options.addOption(flank);
		options.addOption(verbose);
		options.addOption(help);
	}
	
	public static void printHelp() {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp( "AnnotatePairedClusters [options] socrates_paired_cluster_output", options );
	}
	

	private static File tabixFile = null;
	private static MemoryMappedFile tabixMemory = null;
	private static TabixReader annotations = null;

    public static final HashSet<String> SINE_LINE_LTR = new HashSet<String>();
    public static final HashSet<String> SATELLITE = new HashSet<String>();
    public static final HashSet<String> SIMPLE_RPT = new HashSet<String>();
    static {
        SINE_LINE_LTR.add("SINE"); SINE_LINE_LTR.add("LINE"); SINE_LINE_LTR.add("LTR"); SINE_LINE_LTR.add("SINE?"); SINE_LINE_LTR.add("LINE?"); SINE_LINE_LTR.add("LTR?");
        SATELLITE.add("Satellite");
        SIMPLE_RPT.add("Low_complexity"); SIMPLE_RPT.add("Simple_repeat");
    }


	public AnnotatePairedClusters() {
		
	}
	
	public Integer call() {
		return new Integer(0);
	}
	
	public static ArrayList<PairedCluster> loadBreakpoints(String filename) {
		ArrayList<PairedCluster> pb = new ArrayList<PairedCluster>();
		int discarded_chrom = 0;
		try {
			BufferedReader breader = new BufferedReader(new FileReader(filename));
			String line = null; //breader.readLine(); // skip header line
			while ((line=breader.readLine())!=null) {
                if (line.charAt(0)=='#') continue;
				String[] tokens = line.split("\t");
				
				// check chrom
				if (tokens[0].indexOf('_')!=-1 || tokens[3].indexOf('_')!=-1 || tokens[12].indexOf('_')!=-1 || tokens[15].indexOf('_')!=-1) {
					discarded_chrom++;
					continue;
				}
				
                PairedCluster p = new PairedCluster(tokens);
                p.text = line;
				pb.add( p );
			}
			breader.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		System.err.println("Discarded cluster pairs - non-ordinal chromosomes: " + discarded_chrom);
		System.err.println("Loaded cluster pairs: " + pb.size());
		
		return pb;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		CommandLineParser parser = new GnuParser();
		try {
			/* 
			 * Parsing options from command line
			 */
			CommandLine cmd = parser.parse( options, args );
			
            String norm = ((String)cmd.getParsedOptionValue("normal"));
            int flank = cmd.hasOption("flank") ? ((Integer)cmd.getParsedOptionValue("flank")) : 10;
            String rpt = cmd.hasOption("repeatmask") ? ((String)cmd.getParsedOptionValue("repeatmask")) : null;
            System.out.println(norm+"\t"+rpt);
            if (norm==null && rpt==null) {
            	System.err.println("No annotation specified.");
            	System.err.println("Use --normal and/or --repeatmask options to annotate results.");
                printHelp();
                System.exit(1);
            }

            String[] rargs = cmd.getArgs();
			if (rargs.length!=1) {
				printHelp();
				System.exit(1);
			}
           
            ArrayList<PairedCluster> clusterPairs = loadBreakpoints( rargs[0] );
            if (norm != null) annotateNormal( norm, clusterPairs, flank );
            if (rpt != null)  annotateRepeat( rpt, clusterPairs );
            String out = rargs[0]+".annotated";
            outputAnnotation( out, clusterPairs );
		} catch( ParseException exp ) {
	        System.err.println( exp.getMessage() );
	        printHelp();
	        System.exit(1);
	    }
        //ArrayList<PairedCluster> clusterPairs = loadBreakpoints( clustersPairFile );
    }

    public static void annotateNormal(String normalResult, ArrayList<PairedCluster> tumourClusterPairs, int tolerance) {
        TreeSet<PairedCluster> clustersTumour = new TreeSet<PairedCluster>( tumourClusterPairs );
        TreeSet<PairedCluster> clustersNormal = new TreeSet<PairedCluster>( loadBreakpoints( normalResult ) );

		for (PairedCluster tumourCluster : clustersTumour) {
			PairedCluster lower = new PairedCluster(tumourCluster, -2*tolerance);
			PairedCluster upper = new PairedCluster(tumourCluster, 2*tolerance);
			PairedCluster matched = null;
			for (PairedCluster normalCluster : clustersNormal.subSet(lower, true, upper, true)) {
				if (tumourCluster.match(normalCluster, tolerance)) {
					matched = normalCluster;
					break;
				}
			}
			
			if (matched != null) {
				tumourCluster.text += "\tnormal";
			} else {
				tumourCluster.text += "\t";
			}
		}
    }

    public static void annotateRepeat(String rmskFile, ArrayList<PairedCluster> clusterPairs) {
        boolean isBed = rmskFile.indexOf(".bed") != -1;
		
		tabixFile = new File(rmskFile);
		tabixMemory = new MemoryMappedFile(tabixFile, true);
		annotations = new TabixReader(rmskFile, new SeekableMemoryStream(tabixMemory));
		
//		int x=0;
		      
		for (PairedCluster pc : clusterPairs) {
//			x++;
			String ann1 = annotateCluster(pc.cluster1, isBed);
            String ann2 = annotateCluster(pc.cluster2, isBed);

            pc.text += "\t" + ann1 + "\t" + ann2;
		}

	}

    public static void outputAnnotation(String outputFilename, ArrayList<PairedCluster> clusterPairs) {
        try { 
            PrintWriter out = new PrintWriter(new FileOutputStream(new File( outputFilename ) ) ); 
            for (PairedCluster pc : clusterPairs) {
                out.println( pc.text );
            }
            out.close();
        } catch (Exception e) {e.printStackTrace();}
    }

    public static String annotateCluster(Cluster cluster, boolean isBed) {
        String ann1 = "";
        try {
        int s = (cluster.realignPos-10) <= 0 ? 1 : (cluster.realignPos-10);
        String reg = cluster.realignChr+":"+s+"-"+(cluster.realignPos+10);
        TabixReader.Iterator annIter = annotations.query( reg );
        while (annIter!=null && (ann1=annIter.next())!=null) {
            String[] tokens = ann1.split("\t");
            if (!isBed) {
                String rpt = tokens[11];
                if (SINE_LINE_LTR.contains(rpt) || SATELLITE.contains(rpt) || SIMPLE_RPT.contains(rpt)) return rpt;
            } else return tokens[3];
        }} catch (Exception e) { e.printStackTrace(); }
        return "";
    }

    public static String[] toArray(ArrayList<String> list) {
        String[] s = new String[list.size()];
        for (int i=0;i<list.size();i++) { s[i] = list.get(i); }
        return s;
    }

}

class Cluster implements Comparable<Cluster> {
	private String smallChr="", largeChr=""; 
	private int smallPos=-1, largePos=-1, smallDir=-1, largeDir=-1;
	public String realignChr="", anchorChr="";
	public int realignPos=-1, anchorPos=-1;
	public int supportLong=0, supportLongBases, supportShort=0, supportShortBases, supportShortMaxLen;
	public float avg_mapq = 0f;
	
	public boolean realignFwd=false, anchorFwd=false, isShortCluster=false;
	public byte[] realignSeq=null, anchorSeq=null;
	
	public Cluster() {}
	
	public Cluster(String[] fields) {
		String[] realignCoords = fields[0].split(":");
		realignFwd = fields[1].equals("+");
		realignChr = realignCoords[0];
		realignPos = Integer.parseInt(realignCoords[1]);
		String f1 = fields[2].replaceAll("\\*", "");
		realignSeq = Utilities.stringToByte(f1);
		
		String[] anchorCoords = fields[3].split(":");
		anchorFwd = fields[4].equals("+");
		anchorChr = anchorCoords[0];
		anchorPos = Integer.parseInt(anchorCoords[1]);
		String f2 = fields[5].replaceAll("\\*", "");
		anchorSeq = Utilities.stringToByte(f2);

        if (realignSeq.length==0) isShortCluster=true;
		
		this.supportLong = Integer.parseInt(fields[6]);
		this.supportLongBases = Integer.parseInt(fields[7]);
		this.supportShort = Integer.parseInt(fields[8]);
		this.supportShortBases = Integer.parseInt(fields[9]);
		this.supportShortMaxLen = Integer.parseInt(fields[10]);
		this.avg_mapq = Float.parseFloat(fields[11]);

		setSmallLarge();
	}
	
	public void setRealignCoord(String chrom, int pos, boolean fwd) {
		realignChr = chrom; realignPos = pos; realignFwd = fwd; setSmallLarge();
	}
	
	public void setAnchorCoord(String chrom, int pos, boolean fwd) {
		anchorChr = chrom; anchorPos = pos; anchorFwd = fwd; setSmallLarge();
	}
	
	private void setSmallLarge() {
		if (realignChr.compareTo(anchorChr)<0) {
			smallChr = realignChr; smallPos = realignPos; smallDir = realignFwd ? -1 : 1; 
			largeChr = anchorChr; largePos = anchorPos; largeDir = anchorFwd ? -1 : 1;
		} else if (realignChr.compareTo(anchorChr) > 0){
			largeChr = realignChr; largePos = realignPos; largeDir = realignFwd ? -1 : 1;
			smallChr = anchorChr; smallPos = anchorPos; smallDir = anchorFwd ? -1 : 1;
		} else {
			if (realignPos < anchorPos) {
				smallChr = realignChr; smallPos = realignPos; smallDir = realignFwd ? -1 : 1; 
				largeChr = anchorChr; largePos = anchorPos; largeDir = anchorFwd ? -1 : 1;
			} else if (realignPos > anchorPos) {
				largeChr = realignChr; largePos = realignPos; largeDir = realignFwd ? -1 : 1;
				smallChr = anchorChr; smallPos = anchorPos; smallDir = anchorFwd ? -1 : 1;
			} else {
				if (realignFwd && !anchorFwd) {
					smallChr = realignChr; smallPos = realignPos; smallDir = realignFwd ? -1 : 1; 
					largeChr = anchorChr; largePos = anchorPos; largeDir = anchorFwd ? -1 : 1;
				} else {
					largeChr = realignChr; largePos = realignPos; largeDir = realignFwd ? -1 : 1;
					smallChr = anchorChr; smallPos = anchorPos; smallDir = anchorFwd ? -1 : 1;
				}
			}
		}
	}
	
	public String toString(SAMFileInfo info) {
		String c1 = realignChr+":"+realignPos+"\t"+(realignFwd?'+':'-');
		String c2 = anchorChr+":"+anchorPos+"\t"+(anchorFwd?'+':'-');
		return c1 + "\t" + c2 + "\t" + supportLong + "\t" + supportShort;
	}
	
	public int compareTo(Cluster other) {
		if (!smallChr.equals(other.smallChr)) return smallChr.compareTo(other.smallChr);
		if (smallPos != other.smallPos) return smallPos - other.smallPos;
		if (smallDir != other.smallDir) return smallDir - other.smallDir;
		if (!largeChr.equals(other.largeChr)) return largeChr.compareTo(other.largeChr);
		if (largePos != other.largePos) return largePos - other.largePos;
		if (largeDir != other.largeDir) return largeDir - other.largeDir; 
		
		return 0;
	}
	
	public boolean match(Cluster other, int tolerance) {
		if (!smallChr.equals(other.smallChr) || !largeChr.equals(other.largeChr)) return false;
		if (smallDir != other.smallDir || largeDir != other.largeDir) return false;
		return (Math.abs(smallPos - other.smallPos) <= tolerance && Math.abs(largePos - other.largePos) <= tolerance);
	}
}

class PairedCluster implements Comparable<PairedCluster> {
	String text;
	Cluster cluster1, cluster2;
	private Cluster small, large;
	
	
	public PairedCluster(String[] tokens) {
		String[] c1 = Arrays.copyOfRange(tokens, 0, 12);
		String[] c2 = Arrays.copyOfRange(tokens, 12, 24);
		
		cluster1 = new Cluster(c1);
		cluster2 = new Cluster(c2);
		
		//StringBuilder b = new StringBuilder();
		//for (int i=0; i< tokens.length; i++) {
		//	b.append(tokens[i]);
		//	if (i!=tokens.length-1) b.append("\t");
		//}
		//text = b.toString();
		
		if (cluster1.compareTo(cluster2)<0) {
			small = cluster1; large = cluster2;
		} else {
			small = cluster2; large = cluster1;
		}
	}
	
	public PairedCluster(PairedCluster template, int offset) {
		cluster1 = new Cluster();
		cluster1.setRealignCoord(template.cluster1.realignChr, template.cluster1.realignPos+offset, template.cluster1.realignFwd);
		cluster1.setAnchorCoord(template.cluster1.anchorChr, template.cluster1.anchorPos+offset, template.cluster1.anchorFwd); 
		cluster1.realignSeq = template.cluster1.realignSeq;
		cluster1.anchorSeq = template.cluster1.anchorSeq;
		cluster2 = new Cluster();
		cluster2.setRealignCoord(template.cluster2.realignChr, template.cluster2.realignPos+offset, template.cluster2.realignFwd);
		cluster2.setAnchorCoord(template.cluster2.anchorChr, template.cluster2.anchorPos+offset, template.cluster2.anchorFwd); 
		cluster2.realignSeq = template.cluster2.realignSeq;
		cluster2.anchorSeq = template.cluster2.anchorSeq;
        if (cluster2.realignSeq.length==0) cluster2.isShortCluster=true;
		if (cluster1.compareTo(cluster2)<0) {
			small = cluster1; large = cluster2;
		} else {
			small = cluster2; large = cluster1;
		}
	}
	
	public boolean match(PairedCluster other, int tolerance) {
		return (small.match(other.small, tolerance) && large.match(other.large, tolerance));
	}
	
	public int compareTo(PairedCluster other) {
		int smallVsSmall = this.small.compareTo(other.small);
		if (smallVsSmall != 0) return smallVsSmall;
		else return this.large.compareTo(other.large);
	}

    @Override
    public boolean equals(Object o) {
        PairedCluster other = (PairedCluster)o;
        return this.compareTo(other)==0;
    }
	
	public String toString(SAMFileInfo info) {
		return cluster1.toString(info) + "\t" + cluster2.toString(info);
	}
}
