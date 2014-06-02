package au.edu.wehi.idsv;

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
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Rule;
import org.junit.rules.TemporaryFolder;

import au.edu.wehi.idsv.AnnotateBreakends;
import au.edu.wehi.idsv.ClusterEvidence;
import au.edu.wehi.idsv.ExtractEvidence;
import au.edu.wehi.idsv.FileNamingConvention;
import au.edu.wehi.idsv.GenerateDirectedBreakpoints;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.RelevantMetrics;
import au.edu.wehi.idsv.SocratesVariantContext;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;

import com.google.common.collect.Lists;
import com.google.common.io.Files;

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
	public File reference;
	@Before
	public void setup() throws IOException {
		reference = SMALL_FA_FILE;
		input = testFolder.newFile("input.bam");
	}
	public void setReference(File ref) {
		reference = ref;
	}
	public ProcessingContext getCommandlineContext() {
		File metricsFile = new File(testFolder.getRoot(), "input.bam.idsv.working/input.bam.idsv.metrics.txt");
		RelevantMetrics metrics = metricsFile.exists() ? new RelevantMetrics(metricsFile) : null;
		return new ProcessingContext(ReferenceSequenceFileFactory.getReferenceSequenceFile(reference), metrics);
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
	public void createInput(File file) throws IOException {
		if (!file.exists()) throw new IllegalArgumentException(String.format("Missing %s", file));
		Files.copy(file, input);
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
		public int go(boolean perChr) {
			TMP_DIR = Lists.newArrayList(testFolder.getRoot());
			REFERENCE = reference;
			INPUT = input;
			PER_CHR = perChr;
			return doWork();
		}
	}
	public int extractEvidence() {
		return new TestExtractEvidence().go(true);
	}
	public int extractEvidence(boolean perChr) {
		return new TestExtractEvidence().go(perChr);
	}
	public class TestGenerateDirectedBreakpoints extends GenerateDirectedBreakpoints {
		@Override
		protected int doWork() {
			return super.doWork();
		}
		public int go(String chr, Integer minMapq, Integer longsc, Float minId, Float minQual, Integer k) {
			if (minMapq != null) MIN_MAPQ = minMapq;
			if (longsc != null) MIN_BREAKEND_REALIGN_LENGTH = longsc;
			if (minId != null) MIN_PERCENT_IDENTITY = minId;
			if (minQual != null) MIN_LONG_SC_BASE_QUALITY = minQual;
			if (k != null) KMER = k;
			TMP_DIR = Lists.newArrayList(testFolder.getRoot());
			try {
				REFERENCE = reference;
				if (StringUtils.isNotEmpty(chr)) {
					SV_INPUT = FileNamingConvention.getSVBamForChr(input, chr);
					MATE_COORDINATE_INPUT = FileNamingConvention.getMateBamForChr(input, chr);
					VCF_OUTPUT = FileNamingConvention.getBreakendVcfForChr(input, chr);
					fastq = new File(FileNamingConvention.getSVBamForChr(input, chr) + ".fq");
				} else {
					SV_INPUT = FileNamingConvention.getSVBam(input);
					MATE_COORDINATE_INPUT = FileNamingConvention.getMateBam(input);
					VCF_OUTPUT = FileNamingConvention.getBreakendVcf(input);
					fastq = new File(FileNamingConvention.getSVBam(input) + ".fq");
				}
				FASTQ_OUTPUT = fastq;
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			return doWork();
		}
	}
	public int generateDirectedBreakpoints() {
		return new TestGenerateDirectedBreakpoints().go(null, null, null, null, null, null);
	}
	public int generateDirectedBreakpoints(Integer minMapq, Integer longsc, Float minId, Float minQual, Integer k) {
		return new TestGenerateDirectedBreakpoints().go(null, minMapq, longsc, minId, minQual, k);
	}
	public int generateDirectedBreakpoints(String chr) {
		return new TestGenerateDirectedBreakpoints().go(chr, null, null, null, null, null);
	}
	public int generateDirectedBreakpoints(String chr, Integer minMapq, Integer longsc, Float minId, Float minQual, Integer k) {
		return new TestGenerateDirectedBreakpoints().go(chr, minMapq, longsc, minId, minQual, k);
	}	
	public void workingFileShouldExist(String extension) {
		assertTrue(new File(input.getAbsolutePath() + ".idsv.working", input.getName() + extension).exists());
	}
	public class TestClusterEvidence extends ClusterEvidence {
		@Override
		protected int doWork() {
			return super.doWork();
		}
		public int go(List<String> chr) {
			INPUT = input;
			REFERENCE = reference;
			CHROMSOME = chr;
			TMP_DIR = Lists.newArrayList(testFolder.getRoot());
			return doWork();
		}
	}
	public void clusterEvidence() {
		new TestClusterEvidence().go(Lists.<String>newArrayList());	
	}
	public class TestAnnotateBreakends extends AnnotateBreakends {
		@Override
		protected int doWork() {
			return super.doWork();
		}
		public int go() {
			INPUT = input;
			REFERENCE = reference;
			TMP_DIR = Lists.newArrayList(testFolder.getRoot());
			return doWork();
		}
	}
	public void annotateBreakends() {
		new TestAnnotateBreakends().go();
	}
	public List<SAMRecord> getRecords(String extension) {
		File file = new File(input.getAbsolutePath() + ".idsv.working", input.getName() + extension);
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
	public List<SocratesVariantContext> getVcf(String extension) {
		File file = new File(input.getAbsolutePath() + ".idsv.working", input.getName() + extension);
		assertTrue(file.exists());
		VCFFileReader reader = new VCFFileReader(file);
		List<SocratesVariantContext> list = Lists.newArrayList();
		for (VariantContext r : reader) {
			list.add(SocratesVariantContext.create(getCommandlineContext(), r));
		}
		reader.close();
		return list;
	}
	public List<VariantContextDirectedBreakpoint> getVcfBreaks(String extension) {
		File file = new File(input.getAbsolutePath() + ".idsv.working", input.getName() + extension);
		assertTrue(file.exists());
		VCFFileReader reader = new VCFFileReader(file);
		List<VariantContextDirectedBreakpoint> list = Lists.newArrayList();
		for (VariantContext r : reader) {
			SocratesVariantContext vc = SocratesVariantContext.create(getCommandlineContext(), r);
			if (vc instanceof VariantContextDirectedBreakpoint) {
				list.add((VariantContextDirectedBreakpoint)vc);
			}
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
