/**
 * 
 */
package net.wehi.socrates;


import org.apache.commons.lang3.builder.EqualsBuilder;
import org.apache.commons.lang3.builder.HashCodeBuilder;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import java.util.HashSet;
import java.util.HashMap;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.awt.Point;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

import net.sf.samtools.BAMRecordCodec;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordCoordinateComparator;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
//import net.sf.samtools.BAMFileWriter;

import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.SortingCollection;
import net.wehi.socrates.util.GenomicIndexInterval;
import net.wehi.socrates.util.SAMFileInfo;
import net.wehi.socrates.util.SAMRecordSummary;
import net.wehi.socrates.util.RealignmentRecordSummary;
import net.wehi.socrates.util.Utilities;
//import net.wehi.socrates.util.ProgressVerbose;

/**
 * @author hsu
 *
 * Created on Jan 18, 2013
 */

@SuppressWarnings("static-access")
public class BAMStratifier {
	private int scLengthThreshold, scQualityThreshold, minMapq, pcId, bufferSize;
	private float minIdentity;
	private boolean removeDuplicates, disableAllSc; //, removeNonUnique, isBowtie2, isBWA;

	private final int threads;
	private final File sourceBAMFile;
	private final String sourceBAMStem;
	private final SAMFileInfo sourceBAMInfo;
	
	private AtomicInteger counter;
	
	/**
	 * Constructor
	 * Stratify reads into soft-clipped BAM file and anomalous reads BAM file. At the same time, produces gzipped FASTQ files of
	 * long soft-clip part of reads (specified by scLengthThreshold and scQualityThreshold) and one-end-anchored reads.
	 * 
	 * Output filenames have the same stem as source BAM file with suffixes:
	 * _all_sc.bam for soft-clipped alignments
	 * _anomalous.bam for anomalous alignments
	 * _long_sc_l<length>_q<quality>.fastq.gz for long soft-clip sequences
	 * _long_sc_l<length>_q<quality>_anchor.txt.gz for information on unclipped part of soft-clipped alignments
	 * _one_end_anchor.fastq.gz for one-end-anchored sequences.
	 * @param sourceFilename source BAM file
	 * @param number of parallel threads to run
	 * @param scLengthThreshold - minimum length of soft-clipped part
	 * @param scQualityThreshold - minimum average quality of bases of soft-clipped part
	 * @param minMapq - minimum mapping quality of a read to process
	 * @param rmdup - remove duplicated (by alignment position) reads
	 * @param noAllSc - disable outputting of *all_sc.bam files, which are only used for visualization
	 */
	public BAMStratifier(String sourceFilename, int threads, 
						 int scLengthThreshold, int scQualityThreshold, int minMapq, int minPercentIdentity, 
						 boolean rmdup, boolean noAllSc) {
		
		this.threads = threads;
		disableAllSc = noAllSc;
		// template file
		sourceBAMFile = Utilities.assertExists(sourceFilename);
		sourceBAMStem = sourceBAMFile.getName().substring( 0, sourceBAMFile.getName().lastIndexOf('.') );

		// open BAM reader
		SAMFileReader reader = new SAMFileReader( sourceBAMFile );
		reader.setDefaultValidationStringency(SAMFileReader.ValidationStringency.SILENT);
		sourceBAMInfo = new SAMFileInfo(reader);
		// metrics file
		File sourceBAMInfoFile = new File(sourceFilename+".metrics");
		// create metrics file in current directory if one is not found
		if (!sourceBAMInfoFile.exists()) {
			sourceBAMInfoFile = new File(sourceBAMFile.getName()+".metrics");
			if (!sourceBAMInfoFile.exists()) SAMFileInfo.writeBAMMetrics(reader, sourceBAMInfo.libraries, sourceBAMInfoFile.getAbsolutePath(), 10000);

		}
		
		HashMap<String,Integer> readLengths = SAMFileInfo.getLibraryReadLengths(sourceBAMInfoFile.getAbsolutePath());
		int maxReadLength = 0;
		for (Integer length : readLengths.values()) {
			if (length > maxReadLength) maxReadLength = length;
		}
		reader.close();

		this.scLengthThreshold = scLengthThreshold;
		this.scQualityThreshold = scQualityThreshold;
		this.bufferSize = maxReadLength + 10;
		this.removeDuplicates = rmdup;
		this.minIdentity = (float)minPercentIdentity/100f;
		this.pcId = minPercentIdentity;
		this.minMapq = minMapq;
		
		counter = new AtomicInteger(0);
		
		SAMFileWriterFactory.setDefaultCreateIndexWhileWriting(true);
	}
	
	/**
	 * Stratify the entire source BAM file.
	 */
	public void stratifyAll() {
		stratify((ArrayList<GenomicIndexInterval>)null);
	}
	
	/**
	 * Stratify only the given regions. Calls surrogate function stratify(ArrayList<GenomicIndexInterval>).
	 * @param regions Regions in string format. Can contain multiple regions, separated by ";".
	 */
	public void stratify(String regions) {
		if (regions==null) {
			stratify((ArrayList<GenomicIndexInterval>)null);
		} else {
			ArrayList<GenomicIndexInterval> gis = new ArrayList<GenomicIndexInterval>();
			String[] r = regions.split(";", -1);
			for (int i=0; i<r.length; i++) {
				gis.add( GenomicIndexInterval.parse(r[i].trim(), sourceBAMInfo) );
			}
			stratify( gis );
		}
	}
	
	/**
	 * Stratify only the given regions
	 * This is the actual called function for stratify(String regions) and stratifyAll().
	 * All input and output streams are closed just before function returns.
	 * @param regions Regions in ArrayList<GenomicIndexInterval> form.
	 */
	public void stratify(ArrayList<GenomicIndexInterval> regions) {
		try {
			int c = 1;
			ArrayList< Callable<OutputFiles>> tasks = new ArrayList< Callable<OutputFiles>>();
			if (regions==null) {
				for (String s : sourceBAMInfo.sequenceNames) {
					int l = sourceBAMInfo.sequenceLengths.get(s).intValue();
					tasks.add( new StratifyWorker(c, s, 1, l) );
					c++;
				}
			} else {
				for (GenomicIndexInterval gi : regions) {
					tasks.add( new StratifyWorker(c, sourceBAMInfo.getSequenceName(gi.chromIdx), gi.start, gi.end) );
					c++;
				}
			}
			runStratify(tasks);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Invokes multi-threaded execution of tasks.
	 * This is done as a two-step process.
	 * Step 1. stratifies alignments and extract usable SC by reference sequence (StratifyWorker)
	 * Step 2. merge multiple reference sequence BAM files into a single one (MergeWorker)
	 * @param tasks A list of workers to execute in parallel
	 */
	private void runStratify(ArrayList<Callable<OutputFiles>> tasks) {
		ExecutorService pool = null;
		try {
			// start thread service
			pool = Executors.newFixedThreadPool(threads);
			
			// executing concurrent tasks
			java.util.List< Future<OutputFiles> > results = pool.invokeAll(tasks);
			
			int finished = 0;
			int t = tasks.size();
			ArrayList<File> anomalous = new ArrayList<File>(), allSc = new ArrayList<File>(), shortSc = new ArrayList<File>();
			ArrayList<File> longSc = new ArrayList<File>(), oneEndAnchor = new ArrayList<File>();
			File anomalousBAM = new File(sourceBAMStem+"_anomalous.bam"),
				 allScBAM = new File(sourceBAMStem+"_all_sc.bam"),
				 shortScBAM = new File(sourceBAMStem+"_short_sc.bam"),
				 longScFq = new File(sourceBAMStem+"_long_sc_l"+scLengthThreshold+"_q"+scQualityThreshold+"_m"+minMapq+"_i"+pcId+".fastq.gz"),
				 oeaFq = new File(sourceBAMStem+"_one_end_anchor.fastq.gz");
			
			int totAnomalous = 0, totAllSc = 0, totShortSc = 0, totLongSc = 0, totOEA = 0;
			// fetch results from concurrent futures
			for (int i=0; i<t; i++) {
				Future<OutputFiles> future = results.get(i);
				if (future.isDone()) finished++;
				
				OutputFiles outputs = future.get();
				if (outputs.anomalous != null) {
					anomalous.add( outputs.anomalous );
					totAnomalous += outputs.counterAnomalous;
				}
				if (outputs.allSc != null) {
					allSc.add( outputs.allSc );
					totAllSc += outputs.counterAllSc;
				}
				if (outputs.shortSc != null) {
					shortSc.add( outputs.shortSc );
					totShortSc += outputs.counterShortSc;
				}
				if (outputs.longSc != null) {
					longSc.add( outputs.longSc );
					totLongSc += outputs.counterLongSc;
				}
				if (outputs.oneEndAnchor != null) {
					oneEndAnchor.add( outputs.oneEndAnchor );
					totOEA += outputs.counterOneEndAnchor;
				}
			}
			
			if (SOCRATES.verbose) {
				System.err.println(counter.get() + " reads processed");
				System.err.println(totShortSc + " short SC reads");
				System.err.println(totLongSc + " long SC reads");
				System.err.println(totAllSc + " all SC reads");
				System.err.println(totAnomalous + " anomalous reads");
				System.err.println(totOEA + " one-end-anchor reads");
				System.err.println(finished + " straytify tasks finished");
			}
			
			// merge result files
			counter.set(0);
			ArrayList<Callable<File>> mergeTasks = new ArrayList<Callable<File>>();
			mergeTasks.add( new MergeWorker(anomalousBAM, "BAM", anomalous) );
			if (!disableAllSc) mergeTasks.add( new MergeWorker(allScBAM, "BAM", allSc) );
			mergeTasks.add( new MergeWorker(shortScBAM, "BAM", shortSc) );
			mergeTasks.add( new MergeWorker(longScFq, "GZ", longSc));
			mergeTasks.add( new MergeWorker(oeaFq, "GZ", oneEndAnchor));
			
			// executing concurrent tasks
			java.util.List< Future<File> > mergeResults = pool.invokeAll(mergeTasks);
			finished = 0;
			for (int i=0; i<mergeTasks.size(); i++) {
				Future<File> future = mergeResults.get(i);
				if (future.isDone()) finished++;
			}
			if (SOCRATES.verbose) {
				System.err.println(finished + " merge tasks finished");
			}
		} catch (InterruptedException ie) {
			ie.printStackTrace();
		} catch (ExecutionException ee) {
			ee.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			// terminate thread service
			pool.shutdown();
		}
	}
	
	
	class StratifyWorker implements Callable<OutputFiles> {
		public final int threadID;
		private final SAMFileReader sourceBAM;
		private ArrayDeque<HashSet<AlignmentInfo>> hashBuffer = new ArrayDeque<HashSet<AlignmentInfo>>();
		private ArrayDeque<Point> posBuffer = new ArrayDeque<Point>();
		private ArrayDeque<ArrayDeque<SAMRecord> > alnBuffer = new ArrayDeque<ArrayDeque<SAMRecord> >();
		private OutputFiles outputs;
		private final SAMRecordIterator iterator;
		private final String seqName;
		
		public StratifyWorker(int threadID, String seq, int start, int end) throws Exception {
			this.threadID = threadID;
			seqName = seq;
			// create new file reader
			sourceBAM = new SAMFileReader( sourceBAMFile );
			sourceBAM.setDefaultValidationStringency(SAMFileReader.ValidationStringency.SILENT);
			this.iterator = sourceBAM.query(seq, start, end, false);
		}
		
		public OutputFiles call() throws Exception {
			// prepare output files
			outputs = new OutputFiles(seqName, sourceBAM.getFileHeader());
			
			int c = 0;
			while (iterator.hasNext()) {
				SAMRecord aln = iterator.next();
				
				c = counter.incrementAndGet();
				if (c % 50000000 == 0 && SOCRATES.verbose) {
					System.err.println(c + " reads processed");
				}
				
				// skip secondary alignments (important for BWA-MEM aligned files)
				if (aln.getNotPrimaryAlignmentFlag()) continue;
				
				// skip low mapping quality alignments
				if (aln.getMappingQuality() < minMapq || aln.getReadFailsVendorQualityCheckFlag())
					continue;
				
				if (removeDuplicates && aln.getDuplicateReadFlag()) continue;

				// read is NOT mapped, but mate is (One-End-Anchor)
				if (aln.getReadPairedFlag() && aln.getReadUnmappedFlag() && !aln.getMateUnmappedFlag()) {
					writeToFastq(outputs.oneEndAnchoredFASTQ, aln.getReadName()+"&"+aln.getMateReferenceName()+":"+aln.getMateAlignmentStart(),
								 aln.getReadString(), aln.getBaseQualityString());
					outputs.counterOneEndAnchor++;
					continue;
				}
				
				// skip if it's proper and not soft-clipped
				if (aln.getReadPairedFlag()) { 
					if (aln.getProperPairFlag() && !SAMRecordSummary.isAlignmentSoftClipped(aln)) continue;
				} else {
					if (!SAMRecordSummary.isAlignmentSoftClipped(aln)) continue;
				}
				
				addAlignmentToBuffer(aln);
			}
			
			iterator.close();
			close();

			return outputs;
		}
		

		/**
		 * Helper function to write sequence to gzipped FASTQ file.
		 * @param writer
		 * @param header
		 * @param scSeq
		 * @param qual
		 */
		private void writeToFastq(OutputStreamWriter writer, String header, String scSeq, String qual) {
			try {
				writer.write( "@" + header + "\n" );
				writer.write( scSeq + "\n" );
				writer.write( "+\n" );
				writer.write( qual + "\n" );
			} catch (Exception e) {
				e.printStackTrace();
				System.exit(1);
			}
		}
		
		/**
		 * Add an alignment to the buffer. The alignment will be checked against buffered alignments for duplication, only unique fragments are kept.
		 * Buffer auto-flushes and write to BAM file periodically.
		 * @param aln Alignment record to be added to buffer
		 */
		private void addAlignmentToBuffer(SAMRecord aln) {
			AlignmentInfo info = new AlignmentInfo(aln);
			
			// check if a new buffer entry is needed
			if ( this.posBuffer.size() == 0 ) {
				this.hashBuffer.addLast( new HashSet<AlignmentInfo>() );
				this.hashBuffer.getLast().add( info );
				this.alnBuffer.addLast( new ArrayDeque<SAMRecord>() );
				this.alnBuffer.getLast().addLast( aln );
				this.posBuffer.addLast( new Point(aln.getReferenceIndex(), aln.getUnclippedStart()) );
				return;
			}
			
			// new chromosome, flush buffer and add new entry
			if ( this.posBuffer.getLast().x != aln.getReferenceIndex() ) {
				flush();
				this.hashBuffer.addLast( new HashSet<AlignmentInfo>() );
				this.hashBuffer.getLast().add( info );
				this.alnBuffer.addLast( new ArrayDeque<SAMRecord>() );
				this.alnBuffer.getLast().addLast( aln );
				this.posBuffer.addLast( new Point(aln.getReferenceIndex(), aln.getUnclippedStart()) );
				return;
			}
			
			// otherwise, check for duplicate
			if (removeDuplicates) {
				Iterator<HashSet<AlignmentInfo> > it = this.hashBuffer.descendingIterator();
				while (it.hasNext()) {
					HashSet<AlignmentInfo> infoHashSet = it.next();
					if (infoHashSet.contains(info)) return;
				}
			}
		
		    // no duplicate found, append to buffer
			Point lastPos = this.posBuffer.getLast();
			HashSet<AlignmentInfo> lastHash = this.hashBuffer.getLast();
			ArrayDeque<SAMRecord> lastData = this.alnBuffer.getLast();
			if ( lastPos.y == aln.getUnclippedStart() ) {
				lastHash.add( info );
				lastData.add( aln );
			} else {
				this.hashBuffer.addLast( new HashSet<AlignmentInfo>() );
				this.hashBuffer.getLast().add( info );
				this.alnBuffer.addLast( new ArrayDeque<SAMRecord>() );
				this.alnBuffer.getLast().addLast( aln );
				this.posBuffer.add( new Point(aln.getReferenceIndex(), aln.getUnclippedStart()) );
				flush( bufferSize );
			}
			return;
		}


		/**
		 * Flush all alignments in buffer to file.
		 */
		private void flush() {
			flush(0);
		}
		
		/**
		 * Flush alignments to file until only "keepInBuffer" alignments are left in buffer.
		 * @param keepInBuffer
		 */
		private void flush(int keepInBuffer) {
			while (this.posBuffer.size() > keepInBuffer) {
				posBuffer.pollFirst();
				hashBuffer.pollFirst();
				ArrayDeque<SAMRecord> data = alnBuffer.pollFirst();
				
				Iterator<SAMRecord> it = data.iterator();
				while (it.hasNext()) {
					SAMRecord aln = it.next();
					SAMRecordSummary summary = new SAMRecordSummary(aln);

					if (summary.getAlignedPercentIdentity() < minIdentity) continue;
					
					// add to anomalous alignment BAM
					if (aln.getReadPairedFlag() && !(aln.getProperPairFlag())) {
						outputs.anomalousBAM.addAlignment(aln);
						outputs.counterAnomalous++;
					}
					
					if (summary.isClipped()) {
						if (outputs.allScBAM!=null) {
							outputs.allScBAM.addAlignment( aln );
							outputs.counterAllSc++;
						}
					} else { // read is NOT soft clipped at either end
						continue;
					}
					
					boolean proper = aln.getReadPairedFlag() ? aln.getProperPairFlag() : false;
					
					byte[] head = summary.getHeadClipSequence();
					if (head.length >= scLengthThreshold && summary.getAvgHeadClipQuality() >= scQualityThreshold) {
						boolean clip3p = (aln.getReadNegativeStrandFlag()) ? true : false;
						boolean ideal = ((proper&&!clip3p) || (!proper&&clip3p));
						
						String header = RealignmentRecordSummary.makeFastqHeader(aln, summary, aln.getReferenceIndex(), summary.getHeadClipPos(), ideal,
								Utilities.sequenceByteToString(summary.getHeadClipAlignedSequence(), false), '-');
						String clipSeq = Utilities.sequenceByteToString(head, false); // make breakpoint at 3' end
						String clipQual = Utilities.qualityByteToString(summary.getHeadClipQuality(), false);
						writeToFastq(outputs.longScFASTQ, header, clipSeq, clipQual);
						outputs.counterLongSc++;
					}

					byte[] tail = summary.getTailClipSequence();
					if (tail.length >= scLengthThreshold && summary.getAvgTailClipQuality() >= scQualityThreshold) {
						boolean clip3p = (aln.getReadNegativeStrandFlag()) ? false : true;
						boolean ideal = ((proper&&!clip3p) || (!proper&&clip3p));

						String header = RealignmentRecordSummary.makeFastqHeader(aln, summary, aln.getReferenceIndex(), summary.getTailClipPos(), ideal,
								Utilities.sequenceByteToString(summary.getTailClipAlignedSequence(), false), '+');
						String clipSeq = Utilities.sequenceByteToString(tail, false); // make breakpoint at 3' end
						String clipQual = Utilities.qualityByteToString(summary.getTailClipQuality(), false);
						writeToFastq(outputs.longScFASTQ, header, clipSeq, clipQual);
						outputs.counterLongSc++;
					}
					
					if ((head.length>0 && head.length<scLengthThreshold) || (tail.length>0 && tail.length<scLengthThreshold)) {
						outputs.shortScBAM.addAlignment(aln);
						outputs.counterShortSc++;
					}
				}
			}
		}

		/**
		 * Flush all remaining alignments and close all input and output files.
		 */
		private void close() {
			this.flush();
			this.outputs.close();
			this.sourceBAM.close();
		}
	}
	
	class MergeWorker implements Callable<File> {
		File output;
		String format;
		ArrayList<File> sources;
		public MergeWorker(File output, String format, ArrayList<File> sources) {
			this.output = output; this.format = format; this.sources = sources;
		}
		
		public File call() throws Exception {
			if (SOCRATES.verbose) System.err.println("Merging into " + output.getName());
			if (sources.size()==0) {
				if (SOCRATES.verbose) System.err.println("Empty " + format + " " + output.getName());
				return null;
			}

			if (format.equals("BAM")) {
				if (sources.size()==1) {
					sources.get(0).renameTo(output);
					String fn = sources.get(0).getAbsolutePath();
					File indexFile = new File( fn.substring(0,fn.lastIndexOf('.'))+".bai" );
					if (indexFile.exists()) indexFile.renameTo( new File(output.getAbsolutePath()+".bai") );
					indexFile = new File( fn+".bai" );
					if (indexFile.exists()) indexFile.renameTo( new File(output.getAbsolutePath()+".bai") );
				} else {
					SAMFileReader reader = new SAMFileReader( sources.get(0) );
					SAMFileWriterFactory factory = new SAMFileWriterFactory();
					// a sorting buffer
					SortingCollection<SAMRecord> buffer = SortingCollection.newInstance(SAMRecord.class, new BAMRecordCodec(reader.getFileHeader()), 
																						new SAMRecordCoordinateComparator(), 500000, new File("."));
					// output bam file
					SAMFileWriter o = factory.makeBAMWriter( reader.getFileHeader(), true, output );
					reader.close();
					for (File srcFile : sources) {
						SAMFileReader srcReader = new SAMFileReader(srcFile);
						SAMRecordIterator iter = srcReader.iterator();
						while (iter.hasNext()) {
							buffer.add(iter.next());
						}
						iter.close();
						srcReader.close();
						Utilities.deleteBAM(srcFile);
					}
					
					CloseableIterator<SAMRecord> iter2 = buffer.iterator();
					while (iter2.hasNext()) {
						SAMRecord aln = iter2.next();
						o.addAlignment(aln);
					}

					o.close();
				}
			} else /*if (format.equals("GZ"))*/ {
				if (sources.size()==1) {
					sources.get(0).renameTo(output);
				} else {
					if (output.exists()) output.delete();
					PrintWriter o = new PrintWriter(new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream( new FileOutputStream(output) ))));
					for (int i=0; i<sources.size(); i++) {
						File src = sources.get(i);
						BufferedReader r = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(src))));
						String x = null;
						while ((x=r.readLine())!=null) {
							o.println(x);
						}
						r.close();
						sources.get(i).delete();
					}
					o.close();
				}
			}
			
			if (SOCRATES.verbose) System.err.println("Merging into " + output.getName() + " finished");
			return output;
		}
	}
	
	class OutputFiles {
		File anomalous, shortSc, allSc=null, longSc, oneEndAnchor;
		SAMFileWriter allScBAM=null, shortScBAM, anomalousBAM;
		OutputStreamWriter longScFASTQ, oneEndAnchoredFASTQ;
		String seqName;

		int counterAnomalous=0, counterShortSc=0, counterAllSc=0, counterLongSc=0, counterOneEndAnchor=0;
		
		public OutputFiles(String seqName, SAMFileHeader header) throws Exception {
			this.seqName = seqName;
			this.anomalous = new File(sourceBAMStem+"_anomalous_"+seqName+".bam");					//anomalous.deleteOnExit();
			if (!disableAllSc) this.allSc     = new File(sourceBAMStem+"_all_sc_"+seqName+".bam");	//allSc.deleteOnExit();
			this.shortSc   = new File(sourceBAMStem+"_short_sc_"+seqName+".bam");					//shortSc.deleteOnExit();
			this.longSc    = new File(sourceBAMStem+"_long_sc_"+seqName+".fastq.gz");				//longSc.deleteOnExit();
			this.oneEndAnchor = new File(sourceBAMStem+"_one_end_anchor_"+seqName+".fastq.gz");		//oneEndAnchor.deleteOnExit();

			SAMFileWriterFactory factory = new SAMFileWriterFactory();
			// all SC BAM
			if (!disableAllSc) this.allScBAM   = factory.makeBAMWriter( header, false, allSc );
			else this.allScBAM = null;
			// short SC BAM
			this.shortScBAM = factory.makeBAMWriter( header, false, this.shortSc );
			// anomalous BAM
			this.anomalousBAM = factory.makeBAMWriter( header, false, this.anomalous );
			// long SC FASTQ
			this.longScFASTQ = new OutputStreamWriter(new GZIPOutputStream( new FileOutputStream(this.longSc)) );
			// one-end-anchor FASTQ
			this.oneEndAnchoredFASTQ = new OutputStreamWriter(new GZIPOutputStream( new FileOutputStream(this.oneEndAnchor) ));
		}
		
		public void close() {
			try {
				this.anomalousBAM.close();
				this.shortScBAM.close();
				this.longScFASTQ.close();
				this.oneEndAnchoredFASTQ.close();
				if (this.allScBAM != null) this.allScBAM.close();
				delete_empty();
			} catch (Exception e) {
				e.printStackTrace();
				System.exit(1);
			}
		}
		
		/**
		 * Clean up all empty BAM files.
		 */
		public void delete_empty() {
			if (this.counterAnomalous==0) {
				Utilities.deleteBAM(this.anomalous);
				this.anomalous = null;
			}
			if (this.counterShortSc==0) {
				Utilities.deleteBAM(this.shortSc);
				this.shortSc = null;
			}
			if (this.counterAllSc==0 && this.allSc != null) {
				Utilities.deleteBAM(this.allSc);
				this.allSc = null;
			}
			if (this.counterLongSc==0) {
				this.longSc.delete();
				this.longSc = null;
			}
			if (this.counterOneEndAnchor==0) {
				this.oneEndAnchor.delete();
				this.oneEndAnchor = null;
			}
		}
		

	}


	
	private class AlignmentInfo {
		public int refIdx, pos5p, mrefIdx, mpos5p;
		public String strands;
		public boolean singleEnd;
		
		public AlignmentInfo(SAMRecord aln) {
			refIdx = aln.getReferenceIndex();
			pos5p = aln.getReadNegativeStrandFlag() ? aln.getUnclippedEnd() : aln.getUnclippedStart();
		    String strand = aln.getReadNegativeStrandFlag() ? "-" : "+";
		    
			if (aln.getReadPairedFlag()) {
				singleEnd = true;
			    mrefIdx = aln.getMateReferenceIndex();
			    mpos5p = aln.getMateNegativeStrandFlag() ? aln.getMateAlignmentStart()+aln.getReadLength() : aln.getMateAlignmentStart();
			    strands = strand + (aln.getMateNegativeStrandFlag() ? "-" : "+");
			} else {
				singleEnd = false;
				mrefIdx = -1;
				mpos5p = -1;
				strands = strand;
			}
		}
		
		@Override
		public int hashCode() {
			if (this.singleEnd) {
				return new HashCodeBuilder(17,31).
						append( refIdx ).
						append( pos5p ).
						append( strands ).
						toHashCode();
			} else {
				return new HashCodeBuilder(17,31).
							append( refIdx ).
							append( pos5p ).
							append( mrefIdx ).
							append( mpos5p ).
							append( strands ).
							toHashCode();
			}
		}
		
		@Override
		public boolean equals(Object other) {
			if (other==null) return false;
			if (other==this) return true;
			
			AlignmentInfo oinfo = (AlignmentInfo)other;
			if (singleEnd) {
				return new EqualsBuilder().
						append( refIdx, oinfo.refIdx ).
						append( pos5p, oinfo.pos5p ).
						append( strands, oinfo.strands ).
						isEquals();
			} else {
				return new EqualsBuilder().
							append( refIdx, oinfo.refIdx ).
							append( pos5p, oinfo.pos5p ).
							append( mrefIdx, oinfo.mrefIdx ).
							append( mpos5p, oinfo.mpos5p ).
							append( strands, oinfo.strands ).
							isEquals();
			}
		}
	}

    
	private static final Options options = new Options();
	static {
		Option help = new Option( "h", "help", false, "print this message" );
		Option verbose = new Option( "v", "verbose", false, "be verbose of progress" );
		Option noAllSc = new Option( "n", "no-all-sc", false, "do not output *all_sc.bam" );
		Option threads = OptionBuilder.withArgName( "threads" )
										.hasArg()
										.withDescription("Number of threads to use [default: 1]")
										.withType( Number.class )
										.withLongOpt( "threads" )
										.create( 't' );
		Option mapq = OptionBuilder.withArgName( "min-mapq" )
				.hasArg()
				.withDescription("Minimum alignments mapq [default: 5]")
				.withType( Number.class )
				.withLongOpt( "min-mapq" )
				.create( 'q' );
		Option long_sc_len = OptionBuilder.withArgName( "long-sc-len" )
				.hasArg()
				.withDescription("Length threshold of long soft-clip [default: 25 (bp)]")
				.withType( Number.class )
				.withLongOpt( "long-sc-len" )
				.create( 'l' );
		Option keep_duplicate = new Option( "k", "keep-duplicate", false, "keep reads with same alignment position (including mate position for paired-end data");
		Option percent_id = OptionBuilder.withArgName( "percent-id" )
				.hasArg()
				.withDescription("Minimum alignment percent identity to reference [default: 95 (%)]")
				.withType( Number.class )
				.withLongOpt( "percent-id" )
				.create( 'p' );
		Option base_quality = OptionBuilder.withArgName( "base-quality" )
				.hasArg()
				.withDescription("Minimum average base quality score of soft clipped sequence [default: 5]")
				.withType( Number.class )
				.withLongOpt( "base-quality" )
				.create( 'b' );
		
		options.addOption( threads );
		options.addOption( mapq );
		options.addOption( noAllSc );
		options.addOption( long_sc_len );
		options.addOption( percent_id );
		options.addOption( base_quality );
        options.addOption( keep_duplicate );
		options.addOption(verbose);
		options.addOption(help);
	}
	
	public static void printHelp() {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp( "BAMStratifier [options] alignment_bam", options );
	}
	
	public static void main(String[] args) {
		CommandLineParser parser = new GnuParser();
		try {
			/* 
			 * Parsing options from command line
			 */
			CommandLine cmd = parser.parse( options, args );

			int threads = cmd.hasOption("threads") ? ((Number)cmd.getParsedOptionValue("threads")).intValue() : 1;
			int min_mapq = cmd.hasOption("min_mapq") ? ((Number)cmd.getParsedOptionValue("min_mapq")).intValue() : 5;
			int max_long_sc = cmd.hasOption("long-sc-len") ? ((Number)cmd.getParsedOptionValue("long-sc-len")).intValue() : 25;
            int baseq = cmd.hasOption("base-quality") ? ((Number)cmd.getParsedOptionValue("base-quality")).intValue() : 5;
			int min_percent_id = cmd.hasOption("percent-id") ? ((Number)cmd.getParsedOptionValue("percent-id")).intValue() : 95;
			boolean rmdup = cmd.hasOption("keep-duplicate") ? false : true;
			boolean noAllSc = cmd.hasOption("no-all-sc") ? true : false;
            SOCRATES.verbose = cmd.hasOption("verbose");

            String[] remainingArgs = cmd.getArgs();
            if (cmd.hasOption("help") || remainingArgs.length != 1) {
				printHelp();
				System.exit(1);
            }

            String inputBAM = remainingArgs[0];
            Utilities.assertExists(inputBAM);

            System.err.println("\nStratify BAM file " + inputBAM + " with option values:");
            System.err.println("  threads = " + threads);
            System.err.println("  min mapq = " + min_mapq);
            System.err.println("  min long sc length = " + max_long_sc);
            System.err.println("  min avg soft clip base quality = " + baseq);
            System.err.println("  min aligned percent identity = " + min_percent_id);
            System.err.println("  remove duplicates = " + rmdup);
            
       		BAMStratifier b = new BAMStratifier(inputBAM, threads, max_long_sc, baseq, min_mapq, min_percent_id, rmdup, noAllSc);
    	    b.stratifyAll();
        } catch( ParseException exp ) {
	        System.err.println( exp.getMessage() );
	        printHelp();
	        System.exit(1);
	    }

	}
}
