package au.edu.wehi.socrates;

import static org.junit.Assert.assertTrue;
import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.junit.After;
import org.junit.Before;
import org.junit.Rule;
import org.junit.rules.TemporaryFolder;

import com.google.common.collect.Lists;

/**
 * Helper class for testing command-line programs
 * @author Daniel Cameron
 *
 */
public class CommandLineTest extends TestHelper {
	@Rule
    public TemporaryFolder testFolder = new TemporaryFolder();
	public File input;
	public File fastq;
	@Before
	public void setup() throws IOException {
		input = testFolder.newFile("input.bam");
	}
	@After
	public void cleanup() throws Exception {
		input.delete();
		testFolder.delete(); // why isn't @Rule working?
	}
	public void createInput(SAMRecord[] first, SAMRecord[]... data) {
		List<SAMRecord> list = Lists.newArrayList(first);
		for (SAMRecord[] array : data) {
			for (SAMRecord r : array) {
				list.add(r);
			}
		}
		createInput(list.toArray(new SAMRecord[0]));
	}
	public void createInput(SAMRecord... data) {
		SAMFileHeader header = getHeader();
		header.setSortOrder(SortOrder.coordinate);
		SAMFileWriter writer =  new SAMFileWriterFactory().makeSAMWriter(header, true, input);
		SortingCollection<SAMRecord> presort = SortingCollection.newInstance(SAMRecord.class, new BAMRecordCodec(header), new SAMRecordCoordinateComparator(), 100000000, testFolder.getRoot());
		for (SAMRecord r : data) {
			presort.add(r);
		}
		presort.doneAdding();
		for (SAMRecord r : presort) {
			writer.addAlignment(r);
		}
		writer.close();
	}
	public class TestExtractEvidence extends ExtractEvidence {
		@Override
		protected int doWork() {
			return super.doWork();
		}
		public int go() {
			TMP_DIR = Lists.newArrayList(testFolder.getRoot());
			REFERENCE = SMALL_FA_FILE;
			INPUT = input;
			return doWork();
		}
	}
	public int extractEvidence() {
		return new TestExtractEvidence().go();
	}
	public class TestGenerateDirectedBreakpoints extends GenerateDirectedBreakpoints {
		@Override
		protected int doWork() {
			return super.doWork();
		}
		public int go(String chr) {
			return go(chr, null, null, null, null, null);
		}
		public int go(String chr, Integer minMapq, Integer longsc, Float minId, Float minQual, Integer k) {
			if (minMapq != null) MIN_MAPQ = minMapq;
			if (longsc != null) MIN_BREAKEND_REALIGN_LENGTH = longsc;
			if (minId != null) MIN_PERCENT_IDENTITY = minId;
			if (minQual != null) MIN_LONG_SC_BASE_QUALITY = minQual;
			if (k != null) KMER = k;
			TMP_DIR = Lists.newArrayList(testFolder.getRoot());
			try {
				REFERENCE = SMALL_FA_FILE;
				SV_INPUT = FileNamingConvention.getSVBamForChr(input, chr);
				MATE_COORDINATE_INPUT = FileNamingConvention.getMateBamForChr(input, chr);
				VCF_OUTPUT = FileNamingConvention.getBreakendVcfForChr(input, chr);
				fastq = new File(FileNamingConvention.getSVBamForChr(input, chr) + ".fastq");
				FASTQ_OUTPUT = fastq;
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			return doWork();
		}
	}
	public int generateDirectedBreakpoints(String chr) {
		return new TestGenerateDirectedBreakpoints().go(chr);
	}
	public int generateDirectedBreakpoints(String chr, Integer minMapq, Integer longsc, Float minId, Float minQual, Integer k) {
		return new TestGenerateDirectedBreakpoints().go(chr, minMapq, longsc, minId, minQual, k);
	}	
	public void workingFileShouldExist(String extension) {
		assertTrue(new File(input.getAbsolutePath() + ".socrates.working", input.getName() + extension).exists());
	}
	public List<SAMRecord> getRecords(String extension) {
		File file = new File(input.getAbsolutePath() + ".socrates.working", input.getName() + extension);
		assertTrue(file.exists());
		SamReader reader = SamReaderFactory.make().open(file);
		List<SAMRecord> list = Lists.newArrayList();
		for (SAMRecord r : reader) {
			list.add(r);
		}
		try {
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return list;
	}
	public List<VariantContext> getVcf(String extension) {
		File file = new File(input.getAbsolutePath() + ".socrates.working", input.getName() + extension);
		assertTrue(file.exists());
		VCFFileReader reader = new VCFFileReader(file);
		List<VariantContext> list = Lists.newArrayList();
		for (VariantContext r : reader) {
			list.add(r);
		}
		reader.close();
		return list;
	}
	public List<FastqRecord> getFastqRecords() {
		assertTrue(fastq.exists());
		FastqReader reader = new FastqReader(fastq);
		List<FastqRecord> list = Lists.newArrayList();
		for (FastqRecord r : reader) {
			list.add(r);
		}
		reader.close();
		return list;
	}
}
