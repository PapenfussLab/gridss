package au.edu.wehi.socrates;

import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.fastq.FastqReader;
import net.sf.picard.fastq.FastqRecord;
import net.sf.samtools.BAMRecordCodec;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordCoordinateComparator;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.util.SortingCollection;

import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.vcf.VCFFileReader;
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
			if (longsc != null) LONG_SC_LEN = longsc;
			if (minId != null) MIN_PERCENT_IDENTITY = minId;
			if (k != null) MIN_LONG_SC_BASE_QUALITY = minQual;
			if (minQual != null) KMER = k;
			TMP_DIR = Lists.newArrayList(testFolder.getRoot());
			try {
				SV_INPUT = FileNamingConvention.getSVBamForChr(input, chr);
				MATE_COORDINATE_INPUT = FileNamingConvention.getMateBamForChr(input, chr);
				VCF_OUTPUT = FileNamingConvention.getBreakendVcf(input, chr);
				fastq = new File(FileNamingConvention.getSVBamForChr(input, chr) + ".fastq");
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
	public void shouldExist(String file) {
		assertTrue(Files.exists(Paths.get(testFolder.getRoot().toString() + "/" + file)));
	}
	public List<SAMRecord> getRecords(String extension) {
		File file = new File(input.toString() + extension);
		assertTrue(file.exists());
		SAMFileReader reader = new SAMFileReader(file);
		List<SAMRecord> list = Lists.newArrayList();
		for (SAMRecord r : reader) {
			list.add(r);
		}
		reader.close();
		return list;
	}
	public List<VariantContext> getVcf(String extension) {
		File file = new File(input.toString() + extension);
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
