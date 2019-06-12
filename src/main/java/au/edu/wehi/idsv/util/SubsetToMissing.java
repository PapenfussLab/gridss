package au.edu.wehi.idsv.util;

import au.edu.wehi.idsv.ProgressLoggingSAMRecordIterator;
import com.google.common.collect.Lists;
import com.google.common.hash.HashCode;
import com.google.common.hash.HashFunction;
import com.google.common.hash.Hashing;
import htsjdk.samtools.*;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import it.unimi.dsi.fastutil.longs.LongOpenHashBigSet;
import it.unimi.dsi.fastutil.longs.LongSet;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.StandardOptionDefinitions;

import java.io.*;
import java.util.List;
import java.util.Locale;

@CommandLineProgramProperties(
        summary = "Subsets the given BAM file, returning only reads that are not found in the given data set.",  
        oneLineSummary = "Subsets on the given lookup",
        programGroup = gridss.cmdline.programgroups.DataCleaning.class
)
public class SubsetToMissing extends picard.cmdline.CommandLineProgram {
	private Log log = Log.getInstance(SubsetToMissing.class);
	@Argument(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="BAM to test")
    public File INPUT;
	@Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output file")
    public File OUTPUT;
	@Argument(doc="Known reads")
	public List<File> LOOKUP;
	@Argument(doc="Finish processing after this many records", optional=true)
	public Long STOP_AFTER = null;
	@Argument(doc="Size of lookup to preallocate", optional=true)
	public Long PREALLOCATE = null;
	@Override
	protected int doWork() {
		long stop = Long.MAX_VALUE;
		if (STOP_AFTER != null && (long)STOP_AFTER > 0) {
			stop = STOP_AFTER;
		}
		log.debug("Setting language-neutral locale");
    	java.util.Locale.setDefault(Locale.ROOT);
		if (TMP_DIR == null || TMP_DIR.size() == 0) {
			TMP_DIR = Lists.newArrayList(new File("."));
		}
		SamReaderFactory factory = SamReaderFactory.makeDefault()
				.validationStringency(ValidationStringency.SILENT);
		SamReader input = factory.open(INPUT);
		
		
		try (AsyncBufferedIterator<SAMRecord> intputit = new AsyncBufferedIterator<SAMRecord>(input.iterator(), 2, 16384)) {
			
			SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(input.getFileHeader(), true, OUTPUT);
			
			LongSet hashtable;
			if (PREALLOCATE != null) {
				log.info("Preallocating hash table");
				hashtable = new LongOpenHashBigSet(PREALLOCATE);
			} else {
				hashtable = new LongOpenHashBigSet();
			}
			for (File file : LOOKUP) {
				log.info("Loading lookup hashes for " + file.getAbsolutePath());
				SamReader lookup = factory.open(file);
				AsyncBufferedIterator<SAMRecord> it = new AsyncBufferedIterator<SAMRecord>(lookup.iterator(), 2, 16384);
				File cache = new File(file.getPath() + ".SubsetToMissing.cache");
				if (cache.exists()) {
					log.info("Loading lookup hashes from cache");
					long n = stop;
					DataInputStream dis = null;
					try {
						long loadCount = 0;
						dis = new DataInputStream(new BufferedInputStream(new FileInputStream(cache)));
						while (n-- > 0) {
							hashtable.add(dis.readLong());
							if (loadCount % 10000000 == 0) {
								log.info(String.format("Loaded %d from cache", loadCount));
							}
						}
					} catch (EOFException e) {
						try {
							if (dis != null) dis.close();
						} catch (IOException e1) {
							log.error(e1);
						}
					} catch (IOException e) {
						log.error(e);
					}
				} else {
					long n = stop;
					ProgressLoggingSAMRecordIterator loggedit = new ProgressLoggingSAMRecordIterator(it, new ProgressLogger(log));
					try {
						DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(cache)));
						while (loggedit.hasNext() && n-- > 0) {
							long recordhash = hash(loggedit.next());
							hashtable.add(recordhash);
							dos.writeLong(recordhash);
						}
						dos.close();
					} catch (Exception e) {
						log.error(e, "Failed to load lookup. Running with partial results");
					}
					loggedit.close();
				}
				it.close();
			}
			long filtered = 0;
			log.info("Processing input");
			ProgressLoggingSAMRecordIterator it = new ProgressLoggingSAMRecordIterator(intputit, new ProgressLogger(log));
			long n = stop;
			while (it.hasNext() && n-- > 0) {
				SAMRecord r = it.next();
				if (!hashtable.contains(hash(r))) {
					out.addAlignment(r);
				} else {
					filtered++;
					if (filtered % 1000000 == 0) {
						log.info(String.format("Filtered %d reads", filtered));
					}
				}
			}
			log.info("Closing output");
			out.close();
		}
		return 0;
	}
	private HashFunction hf = Hashing.murmur3_128();
	private long hash(SAMRecord r) {
		HashCode hc = hf.newHasher()
			.putBytes(r.getBaseQualities())
			.putBytes(r.getReadBases())
			.hash();
		long hash64 = hc.asLong();
		return hash64;
	}
	public static void main(String[] argv) {
        System.exit(new SubsetToMissing().instanceMain(argv));
    }
}
