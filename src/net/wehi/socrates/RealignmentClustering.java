/**
 * 
 */
package net.wehi.socrates;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.NavigableSet;
import java.util.TreeSet;
import java.util.HashSet;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.awt.geom.Point2D;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.util.SeekableStream;
import net.wehi.socrates.util.GenomicIndexInterval;

//import net.wehi.socrates.util.GenomicIntervalList;
import net.wehi.socrates.util.MemoryMappedFile;
import net.wehi.socrates.util.SAMFileInfo;
import net.wehi.socrates.util.RealignmentRecordSummary;
//import net.wehi.socrates.util.ProgressVerbose;
import net.wehi.socrates.util.RealignmentCluster;
import net.wehi.socrates.util.SAMRecordSummary;
import net.wehi.socrates.util.SeekableMemoryStream;
import net.wehi.socrates.util.Utilities;
import net.sf.picard.sam.BuildBamIndex;

/**
 * @author hsu
 *
 * Created on Jan 24, 2013
 */
@SuppressWarnings("static-access")
public class RealignmentClustering {
	private final String realignScStem;
	private final File realignedScBAMFile, shortScBAMFile;
	private final File realignedScIndexFile, shortScIndexFile;
	private MemoryMappedFile memoryRealignBAM, memoryShortScBAM;
	private MemoryMappedFile memoryRealignIndex, memoryShortScIndex;
	
	private final SAMFileInfo fileInfo;
	private final int minMapq, maxLongSupport, scLenThreshold;
	private final HashMap<String,Float> maxInsertSizes;
	private final int threads, promiscuityThreshold, promiscuityFlank;
	private final float percentIdentityThreshold;
	
	private final boolean threadVerbose = false;
	
	private int jobsTotal=0;
	private AtomicInteger jobsCompleted=new AtomicInteger();
//	private ConcurrentHashMap<Integer,ConcurrentSkipListMap<Integer, AtomicInteger>> promiscuityAnchor;
//	private ConcurrentHashMap<Integer,ConcurrentSkipListMap<Integer, AtomicInteger>> promiscuityRealign;
//	private ConcurrentSkipListSet<>
	
	private static final long MEGABYTE = 1024L * 1024L;

	public static long bytesToMegabytes(long bytes) {
		return bytes / MEGABYTE;
	}

	private boolean idealOnly;
	private boolean find_short_sc_cluster=true;
	
	public RealignmentClustering(String realignedScBAMFilename, 
								 String shortScBAMFilename,
								 String metricsFilename,
								 int threads,
								 int minMapq, int scLenThreshold, int maxSupport, int pidthresh, int pthresh, int pflank, boolean ideal, boolean short_sc_cluster) {
		this.minMapq = minMapq; 
		this.scLenThreshold = scLenThreshold; 
		this.promiscuityThreshold = pthresh; 
		this.promiscuityFlank = pflank;
		this.maxLongSupport = maxSupport; 
		this.percentIdentityThreshold = (float)pidthresh/100f;
		this.idealOnly = ideal;
		this.find_short_sc_cluster = short_sc_cluster;
		this.threads = threads;
		
		realignedScBAMFile = new File(realignedScBAMFilename);
		realignedScIndexFile = findIndexFile(realignedScBAMFilename);
		String fname = realignedScBAMFile.getName();
		realignScStem = fname.substring(0, fname.lastIndexOf('.'));
		
		shortScBAMFile = new File(shortScBAMFilename);
		shortScIndexFile = findIndexFile(shortScBAMFilename);
		
		fileInfo = new SAMFileInfo(realignedScBAMFilename);
		
		HashMap<String,Point2D.Float> insertSizes = SAMFileInfo.getInsertMeanStdev(metricsFilename);
		maxInsertSizes = new HashMap<String, Float>();
		for (String rg : insertSizes.keySet()) {
			Point2D.Float i = insertSizes.get(rg);
			maxInsertSizes.put( rg, i.x + 3*i.y );
		}
		
//		promiscuityAnchor = new ConcurrentHashMap<Integer, ConcurrentSkipListMap<Integer, AtomicInteger>>();
//		promiscuityRealign = new ConcurrentHashMap<Integer, ConcurrentSkipListMap<Integer, AtomicInteger>>();
		
		// load memory mapped BAM file
		memoryRealignBAM = new MemoryMappedFile(realignedScBAMFile, true /* eager loading */);
		memoryRealignIndex = new MemoryMappedFile(realignedScIndexFile, true);
		memoryShortScBAM = new MemoryMappedFile(shortScBAMFile, true /* eager loading */);
		memoryShortScIndex = new MemoryMappedFile(shortScIndexFile, true);		
	}
	
	private File findIndexFile(String bamFilename) {
		File index = new File( bamFilename+".bai" );
		if (!index.exists()) {
			index = new File( bamFilename.substring(0, bamFilename.lastIndexOf("."))+".bai" );
			if (!index.exists()) {
				System.err.println("BAM index file does not exist for " + bamFilename + ".");
				System.err.println("Socrates is creating one for you.");
				BuildBamIndex.createIndex(new SAMFileReader(new File(bamFilename)), index);
			}
		}
		return index;
	}
	
	public void clusterRealignedSC(String regions) throws IOException {
		ArrayList< Callable<ResultCounter> > tasks = new ArrayList<Callable<ResultCounter>>();
		if (regions != null) {
			int id=1;
			String[] r = regions.split(";", -1);
			for (int i=0; i<r.length; i++) {
				GenomicIndexInterval gi1 = GenomicIndexInterval.parse(r[i].trim(), fileInfo);
				for (int j=i; j<r.length; j++) {
					GenomicIndexInterval gi2 = GenomicIndexInterval.parse(r[j].trim(), fileInfo);
					tasks.add( new ClusteringWorker(id, gi1, gi2) );
					id++;
				}
			}
		} else { // no specified region, do clustering chromosome by chromosome
			HashMap<String,Integer> lengths = fileInfo.sequenceLengths;
			
			int id=1;
			for (int i=0; i < fileInfo.sequenceNames.size(); i++) {
				String s = fileInfo.sequenceNames.get(i);
				GenomicIndexInterval gi1 = new GenomicIndexInterval(fileInfo.getSequenceIndex(s), 1, lengths.get(s));
				for (int j=i; j<fileInfo.sequenceNames.size(); j++) {
					String s2 = fileInfo.sequenceNames.get(j);
					GenomicIndexInterval gi2 = new GenomicIndexInterval(fileInfo.getSequenceIndex(s2), 1, lengths.get(s2));
					tasks.add( new ClusteringWorker(id, gi1, gi2) );
					id++;
				}
			}
		}
		
//		progressCounter = new AtomicInteger(0);
		runClustering(tasks);
	}
	
	public void clusterRealignedSC(int blockSize) throws IOException {
		ArrayList< Callable<ResultCounter> > tasks = new ArrayList<Callable<ResultCounter>>();
		ArrayList<GenomicIndexInterval> gis = new ArrayList<GenomicIndexInterval>();
		for (String s : fileInfo.sequenceLengths.keySet()) {
			int sl = fileInfo.sequenceLengths.get(s);
			for (int l=1; l<sl; l+=blockSize) {
				int e = l+blockSize-1;
				if ( e > sl) e = sl; 
				gis.add( new GenomicIndexInterval(fileInfo.getSequenceIndex(s), Math.max(1,l-200), Math.min(sl, e+200)) );
			}
		}
		
		
		int id = 1;
		for (int i=0; i < gis.size(); i++) {
			GenomicIndexInterval gi1 = gis.get(i);
			int sli = fileInfo.getSequenceLength(gi1.chromIdx);
			int lf1 = gi1.start==1 ? 0 : 200, rf1 = gi1.end==sli ? 0 : 200;
			
			for (int j=i; j<gis.size(); j++) {
				GenomicIndexInterval gi2 = gis.get(j);
				int slj = fileInfo.getSequenceLength(gi1.chromIdx);
				int lf2 = gi2.start==1 ? 0 : 200, rf2 = gi2.end==slj ? 0 : 200;
				
				ClusteringWorker worker = new ClusteringWorker(id, gi1, gi2);
				worker.lflank1 = lf1; worker.rflank1 = rf1;
				worker.lflank2 = lf2; worker.rflank2 = rf2;
				tasks.add( worker );
				id++;
			}
		}
			
		runClustering(tasks);
	}
	
	private void runClustering(ArrayList< Callable<ResultCounter> > tasks) {
		ExecutorService pool = null;
		ArrayList<ResultCounter> all_results = new ArrayList<ResultCounter>();
		try {
			// start thread service
			pool = Executors.newFixedThreadPool(threads);

			// executing concurrent tasks
			jobsTotal = tasks.size();
			jobsCompleted.set(0);
			java.util.List< Future<ResultCounter>> results = pool.invokeAll(tasks);
			
			// finishing off
			File pairOut = null;
			pairOut = new File("results_Socrates_paired_"+realignScStem+".txt");
			File unpairOut = new File("results_Socrates_unpaired_"+realignScStem+".txt");
			int pairedLongLong=0, pairedLongShort=0, unpaired=0;
			int finished = 0;

			// merge cluster pairing results
			PrintWriter pairResult = null;
			pairResult = new PrintWriter(pairOut);
			pairResult.println("#C1_realign\tC1_realign_dir\tC1_realign_consensus\tC1_anchor\tC1_anchor_dir\tC1_anchor_consensus\tC1_long_support\tC1_long_support_bases\tC1_short_support\tC1_short_support_bases\tC1_short_support_max_len\tC1_avg_realign_mapq\t" +
							    "C2_realign\tC2_realign_dir\tC2_realign_consensus\tC2_anchor\tC2_anchor_dir\tC2_anchor_consensus\tC2_long_support\tC2_long_support_bases\tC2_short_support\tC2_short_support_bases\tC2_short_support_max_len\tC2_avg_realign_mapq");
			PrintWriter unpairResult = new PrintWriter(unpairOut);
			
			for (int i=0; i<tasks.size(); i++) {
				Future<ResultCounter> future = results.get(i);
				ResultCounter result = future.get();
				all_results.add( result );
				if (future.isDone()) finished++;
				pairedLongLong += result.pairedClusters;
				pairedLongShort += result.longShortPair;
				unpaired += result.unpaired;
				String l = null;
				BufferedReader reader = null;
				if (result.pairResult.exists()) {
					reader = result.getPairResultReader();
					while ((l=reader.readLine())!=null) {
						pairResult.println(l);
					}
					reader.close();
					result.pairResult.delete();
				}
				
				if (result.unpairResult.exists()) {
					reader = result.getUnpairResultReader();
					while ((l=reader.readLine())!=null) {
						unpairResult.println(l);
					}
					reader.close();
					result.unpairResult.delete();
				}
			}
			pairResult.close();
			unpairResult.close();
			
			if (SOCRATES.verbose) {
				System.err.println(finished + " tasks finished");
			}
			System.err.println(pairedLongLong + " clusters paired - long SC to long SC");
			System.err.println(pairedLongShort + " clusters paired - long SC to short SC");
			System.err.println(unpaired + " clusters unpaired");
	    } catch (InterruptedException ie) {
			ie.printStackTrace();
		} catch (ExecutionException ee) {
			ee.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			// terminate thread service
			pool.shutdown();
			// clean up temp files
			for (ResultCounter result : all_results) {
				if (result.pairResult.exists()) result.pairResult.delete();
				if (result.unpairResult.exists()) result.unpairResult.delete();
			}
			Runtime runtime = Runtime.getRuntime();
			long memory = runtime.totalMemory() - runtime.freeMemory();
			System.out.println("Used memory is: "+ bytesToMegabytes(memory) + "MB");
		}
	}
	
	class ClusteringWorker implements Callable<ResultCounter> {
		public int threadID, lflank1=0, rflank1=0, lflank2=0, rflank2=0;
		private GenomicIndexInterval gi1, gi2;
		private SeekableStream realignedScBAMMemoryStream;
		private SeekableStream shortScBAMMemoryStream;
		private SeekableStream realignedScIndexMemoryStream;
		private SeekableStream shortScIndexMemoryStream;
		
		ClusteringWorker(int id, GenomicIndexInterval one, GenomicIndexInterval two) {
			gi1 = one; gi2 = two; threadID = id;
			realignedScBAMMemoryStream = new SeekableMemoryStream(memoryRealignBAM);
			shortScBAMMemoryStream     = new SeekableMemoryStream(memoryShortScBAM);
			realignedScIndexMemoryStream = new SeekableMemoryStream(memoryRealignIndex);
			shortScIndexMemoryStream     = new SeekableMemoryStream(memoryShortScIndex);
		}
		
		public ResultCounter call() throws Exception {
			TreeSet<RealignmentCluster> clusters = new TreeSet<RealignmentCluster>();
			ResultCounter counter = new ResultCounter();
						
			SAMFileReader realignReader = null;
			SAMFileReader shortScReader = null;

			if (threadVerbose) System.err.println("Starting clustering in: " + gi1.toString(fileInfo) + " ---> " + gi2.toString(fileInfo));

			realignReader = new SAMFileReader( realignedScBAMMemoryStream, realignedScIndexMemoryStream, false );
			shortScReader = new SAMFileReader( shortScBAMMemoryStream, shortScIndexMemoryStream, false ); 
			realignReader.setDefaultValidationStringency(SAMFileReader.ValidationStringency.SILENT);
			shortScReader.setDefaultValidationStringency(SAMFileReader.ValidationStringency.SILENT);
			
			// build clusters
			clusters = buildCluster(realignReader, gi1, lflank1, rflank1, gi2, lflank2, rflank2);
			if (!gi1.equals(gi2)) clusters.addAll( buildCluster(realignReader, gi2, lflank2, rflank2, gi1, lflank1, rflank1) );

			// pair clusters
			int pairs = pairLongSCClusters( shortScReader, clusters );
			counter.pairedClusters = pairs;
			
			// pair long SC clusters previously unpaired by looking for short SC cluster
            counter.longShortPair = 0;
            if (find_short_sc_cluster) {
			    pairs = pairLongWithShortClusters( shortScReader, clusters );
			    counter.longShortPair = pairs;
            }
			
			HashSet<RealignmentCluster> paired = new HashSet<RealignmentCluster>();
			int unpair = 0;
			for (RealignmentCluster cluster : clusters) {
				if (cluster.pairedCluster==null) {
					unpair++;
					// output
					counter.unpairResultWriter.println(cluster.toString(fileInfo));
					continue;
				}
				if (cluster.anchorConsensusSeq.length==0) continue;
				if (paired.contains(cluster)) continue;
				paired.add(cluster);
				paired.add(cluster.pairedCluster);
				
				// call breakpoint type
				int c1AnchorToC2Realign = Math.abs( cluster.anchorConsensusPos - cluster.pairedCluster.realignConsensusPos );
				int c2AnchorToC1Realign = Math.abs( cluster.pairedCluster.anchorConsensusPos - cluster.realignConsensusPos );
				String vType = null;
				if (cluster.pairedCluster.anchorConsensusSeq.length==0) vType = "short SC cluster";
				else {
					if (c1AnchorToC2Realign!=c2AnchorToC1Realign) vType = "unequal distances of realigned breakpoint to anchor breakpoint: " + c1AnchorToC2Realign + " v " + c2AnchorToC1Realign; 
					else { // matching distance
						if (c1AnchorToC2Realign==0) { // if distance is zero
							if (cluster.realignConsensusSCSeq.length==0 /*&& cluster.pairedCluster.realignConsensusSCSeq.length==0*/ ) // zero base clipped
								vType = "Blunt-end joining";
							else {
								byte[] c1_realign_sc = cluster.realignConsensusSCSeq;
								byte[] c2_realign_sc = Utilities.getReversedComplementArray(cluster.pairedCluster.realignConsensusSCSeq);
								byte[] longSc = c1_realign_sc.length > c2_realign_sc.length ? c1_realign_sc : c2_realign_sc;
								byte[] shortSc = c1_realign_sc.length > c2_realign_sc.length ? c2_realign_sc : c1_realign_sc;
								if (Utilities.consensusStartsWith(longSc, shortSc, 0.9f)) {
									vType = "Inserted sequence: " + Utilities.sequenceByteToString(longSc, false);
								} else {
									vType = "Unequal inserted sequence: " + Utilities.sequenceByteToString(c1_realign_sc, false) + 
											" v " + Utilities.sequenceByteToString(c2_realign_sc, false);
								}
							}
						} else {
							if (c1AnchorToC2Realign > cluster.anchorConsensusSeq.length || c2AnchorToC1Realign > cluster.pairedCluster.anchorConsensusSeq.length) {
								vType = "Unresolved homology: consensus shorter than length";
							} else {
								byte[] consensus1 = Arrays.copyOfRange(cluster.anchorConsensusSeq, 0, c1AnchorToC2Realign);
								byte[] consensus2 = Arrays.copyOfRange(cluster.pairedCluster.anchorConsensusSeq, 0, c2AnchorToC1Realign);
								Utilities.reverseComplement(consensus2);
								if (Utilities.consensusEqual(consensus1, consensus2, 0.9f)) {
									vType = "Micro-homology: " + c1AnchorToC2Realign + "bp homology found! (" + Utilities.sequenceByteToString(consensus1,false) + ")";
								} else {
									// potentially for equal length imperfect homology!
									vType = "Unresolved sequence difference: " + Utilities.sequenceByteToString(consensus1,false) + 
											" v " + Utilities.sequenceByteToString(consensus2,false);
								}
							}
						}
					}
				}
				
				// output
				String o1 = cluster.toString(fileInfo);
				String o2 = cluster.pairedCluster.toString(fileInfo);
				counter.pairResultWriter.println(o1+"\t"+o2+"\t"+vType);
			}
			counter.unpaired = unpair;
			counter.closeWriters();
			
			realignReader.close();
			shortScReader.close();

			jobsCompleted.incrementAndGet();
            String msg = "";
			if (SOCRATES.verbose) {
                msg += "Thread task: " + gi1.toString(fileInfo) + " ---> " + gi2.toString(fileInfo) + "\n\t";
			}
            msg += "Task completed ( " + jobsCompleted.get() + " / " + jobsTotal + " )";
            System.err.println( msg );
			
			return counter;
		}
		
		private TreeSet<RealignmentCluster> buildCluster(SAMFileReader reader, 
				GenomicIndexInterval realignLocus /*realign locus*/, int realignLeftFlank, int realignRightFlank,
				GenomicIndexInterval anchorLocus /*anchor locus*/, int anchorLeftFlank, int anchorRightFlank) throws Exception {
			TreeSet<RealignmentCluster> clusters = new TreeSet<RealignmentCluster>();

			if (threadVerbose) System.err.println("Building long SC clusters in " + gi1.toString(fileInfo) + " ---> " + gi2.toString(fileInfo));
				
			SAMRecordIterator iter = reader.query(fileInfo.getSequenceName(realignLocus.chromIdx), realignLocus.start, realignLocus.end, false);
			while (iter.hasNext()) {
				SAMRecord aln = iter.next();
				if (aln.getReadUnmappedFlag() || aln.getMappingQuality() < minMapq) continue;
				
				RealignmentRecordSummary summary = new RealignmentRecordSummary(aln);
				
				// skip irrelevant alignments
				if (summary.anchorChr!=anchorLocus.chromIdx) continue;
				if (summary.anchorPos<anchorLocus.start || summary.anchorPos>anchorLocus.end) continue;

				// requires ideal evidence?
				if (idealOnly && !summary.anchorIsIdeal) continue;

				// stringent realignment similarity
				if (summary.calcAlignedPercentIdentity(aln) < percentIdentityThreshold) continue;
				
				// find an existing matching cluster
				RealignmentCluster existing = null;
				RealignmentCluster last = clusters.isEmpty() ? null : clusters.last();
				if (last!=null && last.nearTo(summary, 10)) existing = last;
				else {
					RealignmentCluster lower = new RealignmentCluster(summary, -100);
					RealignmentCluster upper = new RealignmentCluster(summary, 100);
					for (RealignmentCluster cluster : clusters.subSet(lower, true, upper, true)) {
						if (cluster.nearTo(summary, 10) && cluster.matchConsensus(summary, percentIdentityThreshold)) {
							existing = cluster;
							break;
						}
					}
				}
				
				if (existing != null) { // merge if found
					existing.updateWith(summary);
				} else { // add new cluster otherwise
					// add a single-read cluster
					existing = new RealignmentCluster(summary);
					clusters.add( existing );
				}
				
//					progress.stepProgress(prefix, suffix);
			}
			iter.close();
//				progress.end(prefix, suffix);
			
			ArrayList<RealignmentCluster> removable = new ArrayList<RealignmentCluster>();
			// call consensus position
			for (RealignmentCluster cluster : clusters) {
				cluster.callConsensusPosition();
				
				if (cluster.anchorConsensusPos < anchorLocus.start+anchorLeftFlank || cluster.anchorConsensusPos > anchorLocus.end-anchorRightFlank ||
					cluster.realignConsensusPos < realignLocus.start+realignLeftFlank || cluster.realignConsensusPos > realignLocus.end-realignRightFlank) {
					removable.add( cluster );
				}
			}
			
			// remove OB clusters
			for (RealignmentCluster ob : removable) clusters.remove(ob);
			
			removable.clear();
			for (RealignmentCluster cluster : clusters) {
				RealignmentCluster lower = new RealignmentCluster(cluster, -promiscuityFlank), upper = new RealignmentCluster(cluster, promiscuityFlank);
				NavigableSet<RealignmentCluster> sub = clusters.subSet(lower, true, upper, true);
				if (sub.size() > promiscuityThreshold) {
					removable.add( cluster );
					continue;
				}
				cluster.mapqAvg = (float)cluster.mapqTotal / (float)cluster.supportLong;
			}

			// remove promiscuous clusters
			for (RealignmentCluster prom : removable) clusters.remove(prom);

			return clusters;
		}
		
		private int pairLongSCClusters(SAMFileReader reader, TreeSet<RealignmentCluster> clusters) throws Exception {
			if (threadVerbose) System.err.println("Pairing long SC clusters in " + gi1.toString(fileInfo) + " ---> " + gi2.toString(fileInfo));
			
			// find clusters that can be paired
			int p=0;
			HashSet<RealignmentCluster> paired = new HashSet<RealignmentCluster>();
			ArrayList<RealignmentCluster> redundant = new ArrayList<RealignmentCluster>();
			for (RealignmentCluster cluster1 : clusters) {
				if (paired.contains(cluster1)) continue;
				RealignmentCluster lower = new RealignmentCluster(cluster1, -50), upper = new RealignmentCluster(cluster1, 50);
				
				NavigableSet<RealignmentCluster> sub = clusters.subSet(lower, true, upper, true);
				for (RealignmentCluster cluster2 : sub) {
					if (cluster1==cluster2) continue;
					if (paired.contains(cluster2)) {
						redundant.add(cluster1);
						continue;
					}
					if (cluster1.reciprocalNearTo(cluster2, 10)) {
						paired.add(cluster1); paired.add(cluster2);
						cluster1.pairedCluster = cluster2;
						cluster2.pairedCluster = cluster1;
						cluster1.callConsensusSequence();
						cluster2.callConsensusSequence();
						p++;
					}
				}
			}

			//permanently remove redundant clusters from output
//			for (RealignmentCluster redundant_call : redundant) {
//				if(redundant_call.pairedCluster == null) clusters.remove(redundant_call);
//			}

			// add short read support to paired clusters
			for (RealignmentCluster cluster : paired) {
				SAMRecordIterator iter = reader.query(fileInfo.getSequenceName(cluster.anchorChr), 
													  cluster.anchorConsensusPos, cluster.anchorConsensusPos, false);
				
				byte[] realignSeq = Utilities.concatenateByteArrays(cluster.realignConsensusSCSeq, cluster.realignConsensusSeq);
				
				// fetch all short reads
				while (iter.hasNext()) {
					SAMRecord aln = iter.next();
					if (aln.getReadUnmappedFlag() || aln.getMappingQuality() < minMapq) continue;
					
					SAMRecordSummary summary = new SAMRecordSummary(aln);
					byte[] clipSeq = (cluster.anchorForward) ? summary.getTailClipSequence() : summary.getHeadClipSequence();
					
					if (clipSeq.length==0 || clipSeq.length >= scLenThreshold) continue;
					Utilities.reverseComplement(clipSeq); //rev-comp, since clips seq is breakpoint-at-right

					// compare short SC sequence with realigned long SC consensus sequence
					if (Utilities.consensusStartsWith(realignSeq, clipSeq, 0.9f)) {
						cluster.supportShort++;
						cluster.supportShortBases += clipSeq.length;
						if (clipSeq.length > cluster.supportShortMaxLen) cluster.supportShortMaxLen = clipSeq.length;
					}
				}
				iter.close();
			}
			
			return p;
		}
		
		private int pairLongWithShortClusters(SAMFileReader reader, TreeSet<RealignmentCluster> clusters) throws Exception {
			if (threadVerbose) System.err.println("Pairing long SC clusters with short SC clusters in " + gi1.toString(fileInfo) + " ---> " + gi2.toString(fileInfo));
			
			int p=0;
			TreeSet<RealignmentCluster> shortScClusters = new TreeSet<RealignmentCluster>();

			// look for short read support in unpaired clusters
			for (RealignmentCluster cluster : clusters) {
				if (cluster.pairedCluster != null) {
					continue; // skip long-long paired clusters
				}
				if (cluster.supportLong > maxLongSupport) {
					continue; // skip deeply supported clusters where we expect existing reciprocal support
				}
				
				cluster.callConsensusSequence();
				
				int maxScLen = 0;
				int cumulativeScLen = 0;
				int shortReads = 0;
				int cumulativeMapq = 0;
				
				byte[] consensusSeq = cluster.realignConsensusSCSeq.length==0 ? 
						cluster.anchorConsensusSeq : Utilities.getReversedComplementArray(cluster.realignConsensusSCSeq);
				
				SAMRecordIterator iter = reader.query(fileInfo.getSequenceName(cluster.realignChr), 
						  cluster.realignConsensusPos-10, cluster.realignConsensusPos+10, false);
				
				// fetch all short reads
				while (iter.hasNext()) {
					SAMRecord aln = iter.next();
					if (aln.getReadUnmappedFlag() || aln.getMappingQuality() < minMapq) continue;
					
					SAMRecordSummary summary = new SAMRecordSummary(aln);
					int offset = (cluster.realignForward) ? summary.getTailClipPos() - cluster.realignConsensusPos : cluster.realignConsensusPos - summary.getHeadClipPos();
					if (Math.abs(offset)>10) continue;
					
					byte[] clipSeq = (cluster.realignForward) ? summary.getTailClipSequence() : summary.getHeadClipSequence();
					if (clipSeq.length==0 || clipSeq.length >= scLenThreshold) continue;
					Utilities.reverseComplement(clipSeq);
					
					
					byte[] clipSeqFinal=null, consensusFinal=null;
					if (offset > 0) {
						clipSeqFinal = clipSeq;
						if (offset>=consensusSeq.length) consensusFinal=new byte[0];
						else consensusFinal = Arrays.copyOfRange(consensusSeq, offset, consensusSeq.length);
					} else if (offset < 0 ) {
						if (-offset>=clipSeq.length) clipSeqFinal = new byte[0];
						else clipSeqFinal = Arrays.copyOfRange(clipSeq, -offset, clipSeq.length);
						consensusFinal = consensusSeq;
					} else {
						clipSeqFinal = clipSeq;
						consensusFinal = consensusSeq;
					}
					
					if (Utilities.consensusStartsWith(consensusFinal, clipSeqFinal, 0.9f)) {
						if (clipSeqFinal.length > maxScLen) maxScLen = clipSeqFinal.length;
						cumulativeScLen += clipSeqFinal.length;
						cumulativeMapq += aln.getMappingQuality();
						shortReads++;
					}
					
//					if (cumulativeScLen >= 10 && maxScLen >= 5 && shortReads >= 2) {
//						break;
//					}
				}
				iter.close();
				
				if (cumulativeScLen >= 10 && maxScLen >= 5 && shortReads >= 2) {
					// create new cluster and pair it
					RealignmentCluster shortScCluster = new RealignmentCluster(cluster.anchorChr, cluster.anchorConsensusPos, cluster.anchorForward, 
																			   cluster.realignChr, cluster.realignConsensusPos, cluster.realignForward);
					shortScCluster.supportShort = shortReads;
					shortScCluster.supportShortBases = cumulativeScLen;
					shortScCluster.supportShortMaxLen = maxScLen;
					shortScCluster.mapqAvg = (float)cumulativeMapq / (float)shortReads;
					shortScClusters.add(shortScCluster);
					cluster.pairedCluster = shortScCluster;
					
					// add short read support to paired clusters
					iter = reader.query(fileInfo.getSequenceName(cluster.anchorChr), 
														  cluster.anchorConsensusPos, cluster.anchorConsensusPos, false);
						
					byte[] realignSeq = Utilities.concatenateByteArrays(cluster.realignConsensusSCSeq, cluster.realignConsensusSeq);
					
					// fetch all short reads
					while (iter.hasNext()) {
						SAMRecord aln = iter.next();
						if (aln.getReadUnmappedFlag() || aln.getMappingQuality() < minMapq) continue;
						
						SAMRecordSummary summary = new SAMRecordSummary(aln);
						byte[] clipSeq = (cluster.anchorForward) ? summary.getTailClipSequence() : summary.getHeadClipSequence();
						if (clipSeq.length==0 || clipSeq.length >= scLenThreshold) continue;
						Utilities.reverseComplement(clipSeq); //rev-comp, since clips seq is breakpoint-at-right

						// compare short SC sequence with realigned long SC consensus sequence
						if (Utilities.consensusStartsWith(realignSeq, clipSeq, 0.9f)) {
							cluster.supportShort++;
							cluster.supportShortBases += clipSeq.length;
							if (clipSeq.length > cluster.supportShortMaxLen) cluster.supportShortMaxLen = clipSeq.length;
						}
					}
					iter.close();
					
					p++;
				}
			}
			
			clusters.addAll(shortScClusters);
			
			return p;
		}
	}
	
	class ResultCounter {
		int pairedClusters = 0, longShortPair = 0, unpaired = 0;
		PrintWriter pairResultWriter=null;	File pairResult=null;
		PrintWriter unpairResultWriter; File unpairResult;
		public ResultCounter() throws Exception {
			File dir = new File(System.getProperty("user.dir"));
			pairResult = File.createTempFile("Socrates_paired_"+realignScStem, ".tmp", dir);
			pairResult.deleteOnExit();
			pairResultWriter = new PrintWriter(pairResult);
			unpairResult = File.createTempFile("Socrates_unpaired_"+realignScStem, ".tmp", dir);
			unpairResult.deleteOnExit();
			unpairResultWriter = new PrintWriter(unpairResult);
		}
		
		public void closeWriters() {
			pairResultWriter.close();
			// delete empty files
			if (pairedClusters+longShortPair==0) pairResult.delete();
			unpairResultWriter.close();
			// delete empty files
			if (unpaired==0) unpairResult.delete();
		}
		
		public BufferedReader getPairResultReader() throws Exception {
			return new BufferedReader(new FileReader(pairResult));
		}
		
		public BufferedReader getUnpairResultReader() throws Exception {
			return new BufferedReader(new FileReader(unpairResult));
		}
	}
	
	
	private static final Options options = new Options();
	static {
		Option help = new Option( "h", "help", false, "print this message" );
		Option verbose = new Option( "v", "verbose", false, "be very verbose of progress" );
		Option threads = OptionBuilder.withArgName( "threads" )
										.hasArg()
										.withDescription("Number of threads to use [default: 1]")
										.withType( Number.class )
										.withLongOpt( "threads" )
										.create( 't' );
		Option ideal = new Option( "i", "ideal-only", false, "FOR PAIRED-END DATA ONLY. Use only proper pair 5' SC and anomalous pair 3' SC [default: false]");
		Option mapq = OptionBuilder.withArgName( "min-mapq" )
				.hasArg()
				.withDescription("Minimum realignments mapq [default: 2]")
				.withType( Number.class )
				.withLongOpt( "min-mapq" )
				.create( 'q' );
		Option long_sc_len = OptionBuilder.withArgName( "long-sc-len" )
				.hasArg()
				.withDescription("Length threshold of long soft-clip [default: 25 (bp)]")
				.withType( Number.class )
				.withLongOpt( "long-sc-len" )
				.create( 'l' );

		Option max_support = OptionBuilder.withArgName( "max-support" )
				.hasArg()
				.withDescription("Maximum realignment support to search for short SC cluster [default: 30]. --no-short-sc-cluster option cannot be used.")
				.withType( Number.class )
				.withLongOpt( "max-support" )
				.create( 's' );
		Option no_short_sc_cluster = new Option( "c", "no-short-sc-cluster", false, "Disable search for short soft clip cluster support for unpaired clusters.");

		Option percent_id = OptionBuilder.withArgName( "percent-id" )
				.hasArg()
				.withDescription("Minimum realignment percent identity to reference [default: 95 (%)]")
				.withType( Number.class )
				.withLongOpt( "percent-id" )
				.create( 'p' );

		Option promiscuity = OptionBuilder.withArgName( "promiscuity" )
				.hasArg()
				.withDescription("Exclude cluster if more than [promiscuity] clusters within [flank]bp of a breakpoint [default: 5]")
				.withType( Number.class )
				.withLongOpt( "promiscuity" )
				.create( 'm' );
		Option flank = OptionBuilder.withArgName( "flank" )
				.hasArg()
				.withDescription("Size of flank for promiscuity filter [default: 50 (bp)]")
				.withType( Number.class )
				.withLongOpt( "flank" )
				.create( 'f' );
		
		options.addOption( threads );
		options.addOption( ideal );
		options.addOption( mapq );
		options.addOption( long_sc_len );
		options.addOption( max_support );
		options.addOption( percent_id );
		options.addOption( promiscuity );
		options.addOption( flank );
		options.addOption( no_short_sc_cluster );
		options.addOption(verbose);
		options.addOption(help);
	}
	
	public static void printHelp() {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp( "RealignmentClustering [options] realigned_sc_bam short_sc_bam metrics_file [block_size]", options );
	}
	
	public static void main(String[] args) throws IOException {
		CommandLineParser parser = new GnuParser();
		try {
			/* 
			 * Parsing options from command line
			 */
			CommandLine cmd = parser.parse( options, args );
			
			int threads = cmd.hasOption("threads") ? ((Number)cmd.getParsedOptionValue("threads")).intValue() : 1;
			boolean ideal = cmd.hasOption("ideal-only") ? true : false;
			int min_mapq = cmd.hasOption("min_mapq") ? ((Number)cmd.getParsedOptionValue("min_mapq")).intValue() : 2;
			int max_long_sc = cmd.hasOption("long-sc-len") ? ((Number)cmd.getParsedOptionValue("long-sc-len")).intValue() : 25;
			int max_support = cmd.hasOption("max_support") ? ((Number)cmd.getParsedOptionValue("max_support")).intValue() : 30;
			boolean short_sc_cluster = cmd.hasOption("no-short-sc-cluster") ? false : true;
			int min_percent_id = cmd.hasOption("percent-id") ? ((Number)cmd.getParsedOptionValue("percent-id")).intValue() : 95;
			int promiscuity = cmd.hasOption("promiscuity") ? ((Number)cmd.getParsedOptionValue("promiscuity")).intValue() : 5;
			int pflank = cmd.hasOption("flank") ? ((Number)cmd.getParsedOptionValue("flank")).intValue() : 50;
			SOCRATES.verbose = cmd.hasOption("verbose");
			
			String[] remainingArgs = cmd.getArgs();
			if (remainingArgs.length!=3 && remainingArgs.length!=4) {
				printHelp();
				System.exit(1);
			}
			String realignedSc = remainingArgs[0];
			String shortSc = remainingArgs[1];
			String metrics = remainingArgs[2];
			String region = remainingArgs.length>=4 ? remainingArgs[3] : null;

            System.err.println( "\nClustering re-alignments." );
            System.err.println( "  search for short soft clip cluster: " + short_sc_cluster );
            if (short_sc_cluster)
            	System.err.println( "  maximum long re-alignment support for unpaired cluster to search for short soft clip cluster: " + max_support);
            System.err.println( "  ideal evidence only: " + ideal );
            System.err.println( "  parallel threads: " + threads);
            System.err.println( "  minimum mapq of re-alignment: " + min_mapq );
            System.err.println( "  minimum \"long\" soft clip length: " + max_long_sc );
            System.err.println( "  minimum percent identity of realignment: " + min_percent_id );
            System.err.println( "  minimum number of clusters in a region to be promiscuous region: " + promiscuity );
            System.err.println( "  size of region: " + pflank +"bp" );
			
			RealignmentClustering r = new RealignmentClustering(
					realignedSc, 
					shortSc, 
					metrics, 
					threads /*threads*/, 
					min_mapq /*minMapq*/, 
					max_support /*maxSupport*/, 
					max_long_sc,
					min_percent_id,
					promiscuity,
					pflank,
					ideal /*idealOnly*/,
					short_sc_cluster /*use short_sc_cluster*/);

			if (region!=null && region.length()>=1 && region.matches("[0-9]+")) {
				int rsize = Integer.parseInt(region);
                System.err.println( "  using bucket size: " + rsize );
				if (rsize>0) r.clusterRealignedSC(Integer.parseInt(region));
				else r.clusterRealignedSC(null);
			}
			else {
                r.clusterRealignedSC(null);
            }
		} catch( ParseException exp ) {
	        System.err.println( exp.getMessage() );
	        printHelp();
	        System.exit(1);
	    }
	}
}
