package gridss;

import au.edu.wehi.idsv.*;
import au.edu.wehi.idsv.bed.BedpeIterator;
import au.edu.wehi.idsv.bed.BedpeRecord;
import au.edu.wehi.idsv.bed.BedpeWriter;
import au.edu.wehi.idsv.configuration.GridssConfiguration;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import au.edu.wehi.idsv.util.WindowedSortingIterator;
import au.edu.wehi.idsv.vcf.VcfFormatAttributes;
import au.edu.wehi.idsv.vcf.VcfSvConstants;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterators;
import com.google.common.collect.Ordering;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.math3.util.Pair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;


@CommandLineProgramProperties(
		summary = "Generates a PON from multiple input files",
		oneLineSummary = "Generates a PON from multiple input files",
        programGroup = gridss.cmdline.programgroups.DataConversion.class)
public class GeneratePonBedpe extends CommandLineProgram {
	private static final Log log = Log.getInstance(GeneratePonBedpe.class);
    @Argument(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input GRIDSS VCFs", optional=false)
    public List<File> INPUT;
    @Argument(doc="Existing GRIDSS breakpoint PON", optional=true)
	public File INPUT_BEDPE = null;
	@Argument(doc="Existing GRIDSS single breakend PON", optional=true)
	public File INPUT_BED = null;
    @Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output BEDPE", optional=false)
    public File OUTPUT_BEDPE;
	@Argument(shortName="SBO", doc="Output BED", optional=false)
	public File OUTPUT_BED;
	@Argument(shortName="NO", doc="0-based ordinals of the normal samples in the VCF.", optional=false)
    public List<Integer> NORMAL_ORDINAL = new ArrayList<>();
	@Argument(shortName="Q", doc="Minimum variant quality score for a breakpoint variant to be considered part of the normal.", optional=true)
	public double MIN_BREAKPOINT_QUAL = 75;
	@Argument(shortName="BQ", doc="Minimum variant quality score for a single breakend variant to be considered part of the normal.", optional=true)
	public double MIN_BREAKEND_QUAL = 428;
	@Argument(doc="Include imprecise calls in the panel of normals.", optional=true)
	public boolean INCLUDE_IMPRECISE_CALLS = false;
	private int MAX_BREAKPOINT_HOMOLOGY_LENGTH = 2000;
	@Argument(doc="Number of worker threads to spawn. Defaults to number of cores available with a maximum of one 1 per input file.", shortName="THREADS")
	public int WORKER_THREADS = Runtime.getRuntime().availableProcessors();

	@Override
	protected String[] customCommandLineValidation() {
		if (INPUT_BED != null ^ INPUT_BEDPE != null) {
			return new String[] {"INPUT_BED and INPUT_BEDPE must be specified together."};
		}
		if (OUTPUT_BED != null && OUTPUT_BED.equals(INPUT_BED)) {
			return new String[] {"INPUT_BED and OUTPUT_BED must be different files."};
		}
		if (OUTPUT_BEDPE != null && OUTPUT_BEDPE.equals(INPUT_BEDPE)) {
			return new String[] {"INPUT_BEDPE and OUTPUT_BEDPE must be different files."};
		}
		if (NORMAL_ORDINAL == null || NORMAL_ORDINAL.size() == 0) {
			return new String[] {"NORMAL_ORDINAL must be specified"};
		}
		return super.customCommandLineValidation();
	}

	@Override
	protected boolean requiresReference() {
		return true;
	}

    @Override
	protected int doWork() {
		try {
			if (INPUT_BED != null) {
				IOUtil.assertFileIsReadable(INPUT_BED);
			}
			if (INPUT_BEDPE != null) {
				IOUtil.assertFileIsReadable(INPUT_BEDPE);
			}
			IOUtil.assertFileIsWritable(OUTPUT_BEDPE);
			IOUtil.assertFileIsWritable(OUTPUT_BED);
			for (File f : INPUT) {
				IOUtil.assertFileIsReadable(f);
			}
			if (INCLUDE_IMPRECISE_CALLS) {
				log.error("Imprecise call inclusion not recommended due to overly aggressive PON matching.");
				Thread.sleep(10000);
			}
			log.debug("Setting language-neutral locale");
			Locale.setDefault(Locale.ROOT);
			GridssConfiguration config;
			try {
				config = new GridssConfiguration((File)null, TMP_DIR.get(0));
			} catch (ConfigurationException e) {
				throw new RuntimeException(e);
			}
			GenomicProcessingContext pc = new GenomicProcessingContext(new FileSystemContext(TMP_DIR.get(0), TMP_DIR.get(0), MAX_RECORDS_IN_RAM), REFERENCE_SEQUENCE, null);
			pc.setCommandLineProgram(this);
			Iterator<Pair<BreakendSummary, Integer>> mergedIt = Iterators.mergeSorted(ImmutableList.of(
					filteredMerge(pc, INPUT, NORMAL_ORDINAL),
					getExistingPON(pc.getDictionary(), INPUT_BEDPE, INPUT_BED)), ByBreakendStartEnd);
			BedpeMergingCounter pe = new BedpeMergingCounter();
			BedMergingCounter se = new BedMergingCounter(true);
			BedpeWriter writer = new BedpeWriter(pc.getDictionary(), OUTPUT_BEDPE);
			BufferedWriter seWriter = Files.newBufferedWriter(OUTPUT_BED.toPath(), StandardCharsets.US_ASCII);
			while(mergedIt.hasNext()) {
				Pair<BreakendSummary, Integer> record = mergedIt.next();
				if (record.getFirst() instanceof BreakpointSummary) {
					Pair<BreakpointSummary, Integer> bpRecord = Pair.create((BreakpointSummary)record.getFirst(), record.getSecond());
					writeBedpe(pe.process(bpRecord), writer);
				} else {
					writeBed(pc.getReference().getSequenceDictionary(), seWriter, se.process(record));
				}
			}
			writeBedpe(pe.finish(), writer);
			writeBed(pc.getReference().getSequenceDictionary(), seWriter, se.finish());
			writer.close();
			seWriter.close();
		} catch (IOException e) {
			log.error(e);
			return 1;
		} catch (InterruptedException e) {
		}
		log.error("Imprecise call inclusion not recommended due to overly aggressive PON matching.");
		return 0;
 	}
 	private static Pair<BreakendSummary, Integer> toPair(SAMSequenceDictionary dictionary, BEDFeature feat) {
		String chr = feat.getContig();
		int start = feat.getStart();
		int end = feat.getEnd();
		int referenceIndex = dictionary.getSequenceIndex(chr);
		Strand strand = feat.getStrand();
		float score = feat.getScore();
		return Pair.create(new BreakendSummary(referenceIndex, strand == Strand.FORWARD ? BreakendDirection.Forward : BreakendDirection.Backward, start, start, end), (int)score);
	}

	private Iterator<Pair<BreakendSummary, Integer>> getExistingPON(SAMSequenceDictionary dictionary, File bedpeFile, File bedFile) throws IOException {
		if (bedpeFile == null || bedFile == null) {
			return ImmutableList.<Pair<BreakendSummary, Integer>>of().iterator();
		}
		BedpeIterator peIt = new BedpeIterator(bedpeFile, dictionary);
		Iterator<Pair<BreakendSummary, Integer>> peItTransformed = Iterators.transform(peIt, (BedpeRecord x) -> Pair.create(x.bp, Integer.parseInt(x.score)));
		BEDCodec codec = new BEDCodec();
		AbstractFeatureReader<BEDFeature, LineIterator> reader = AbstractFeatureReader.getFeatureReader(bedFile.getPath(), codec, false);
		Iterator<Pair<BreakendSummary, Integer>> seItTransformed = Iterators.transform(reader.iterator(), (BEDFeature x) -> toPair(dictionary, x));
		return new AutoClosingIterator<>(Iterators.mergeSorted(ImmutableList.of(peItTransformed, seItTransformed), ByBreakendStartEnd), reader, peIt);
	}
	private static <T> Iterable<T> expand(Pair<T, Integer> pair) {
		List<T> result = new ArrayList<>(pair.getSecond());
		for (int i = 0; i < pair.getSecond(); i++) {
			result.add(pair.getFirst());
		}
		return result;
	}

	private void writeBed(SAMSequenceDictionary dict, BufferedWriter writer, List<Pair<BreakendSummary, Integer>> list) throws IOException {
		for (Pair<BreakendSummary, Integer> pair : list) {
			writeBed(dict, writer, pair.getFirst(), pair.getSecond());
		}
	}
	private void writeBed(SAMSequenceDictionary dict, BufferedWriter writer, BreakendSummary bs, int count) throws IOException {
		int referenceIndex = bs.referenceIndex;
		int bedStart = bs.start - 1;
		int bedEnd = bs.end;
		writer.write(String.format("%s\t%d\t%d\t.\t%d\t%s\n", dict.getSequence(referenceIndex).getSequenceName(), bedStart, bedEnd, count, bs.direction == BreakendDirection.Forward ? '+' : '-'));
	}

	private void writeBedpe(List<Pair<BreakpointSummary, Integer>> list, BedpeWriter writer) throws IOException {
		for (Pair<BreakpointSummary, Integer> emitted : list) {
			writer.write(emitted.getFirst(), ".", Integer.toString(emitted.getSecond()));
		}
	}

	private Iterator<Pair<BreakendSummary, Integer>> filteredMerge(GenomicProcessingContext pc, List<File> file, List<Integer> ordinals) {
    	List<CloseableIterator<Pair<BreakendSummary, Integer>>> fileIt = new ArrayList<>();
    	for (File f : file) {
    		fileIt.add(getFilteredIterator(pc, f, ordinals));
		}
		// allocate VCF files to worker threads
		final AtomicInteger counter = new AtomicInteger(0);
		final int size = WORKER_THREADS;
		final List<Iterator<Pair<BreakendSummary, Integer>>> partitioned = fileIt.stream()
				.collect(Collectors.groupingBy(it -> counter.getAndIncrement() % size))
				.values()
				.stream()
				.map(list -> list.size() == 1 ? list.get(0) : Iterators.mergeSorted(list, ByBreakendStartEnd))
				.map(it -> new AsyncBufferedIterator<>(it,"AsyncVCF"))
				.collect(Collectors.toList());
		Iterator<Pair<BreakendSummary, Integer>> mergedIt = new AsyncBufferedIterator<>(Iterators.mergeSorted(partitioned, ByBreakendStartEnd), "Merged VCF reader");
		return mergedIt;
	}
	private CloseableIterator<Pair<BreakendSummary, Integer>> getFilteredIterator(GenomicProcessingContext pc, File file, List<Integer> ordinals) {
		VCFFileReader vcfReader = new VCFFileReader(file, false);
		CloseableIterator<VariantContext> it = vcfReader.iterator();
		Iterator<Pair<BreakendSummary, Integer>> idsvIt = Iterators.transform(it, variant -> getBreakendSummary(pc, variant, ordinals));
		Iterator<Pair<BreakendSummary, Integer>> idsvFilteredIt = Iterators.filter(idsvIt, pair -> pair.getSecond() > 0);
		Iterator<Pair<BreakendSummary, Integer>> bpit = new PairBreakendSummaryWindowedSortingIterator(pc.getLinear(), idsvFilteredIt, MAX_BREAKPOINT_HOMOLOGY_LENGTH);
		return new AutoClosingIterator<>(bpit, it, vcfReader);
	}
	private Pair<BreakendSummary, Integer> getBreakendSummary(GenomicProcessingContext pc, VariantContext variant, List<Integer> ordinals) {
		IdsvVariantContext vc = IdsvVariantContext.create(pc, null, variant);
		if (!INCLUDE_IMPRECISE_CALLS && vc.hasAttribute(VcfSvConstants.IMPRECISE_KEY)) {
			return null;
		}
		if (!(vc instanceof DirectedEvidence)) {
			return null;
		}
		int normalCount = passesNormalFilter(vc, ordinals);
		return Pair.create(((DirectedEvidence)vc).getBreakendSummary(), normalCount);
	}

	private int passesNormalFilter(IdsvVariantContext vc, List<Integer> ordinals) {
		int count = 0;
		for (int ordinal : ordinals) {
			if (passesNormalFilter(vc, ordinal)) {
				count++;
			}
		}
		return count;
	}
	private boolean passesNormalFilter(IdsvVariantContext vc, int ordinal) {
		Genotype geno = vc.getGenotype(ordinal);
		if (geno == null) {
			throw new RuntimeException(String.format("Missing ordinal %d for ", ordinal, vc.getID()));
		}
		double bpqual = Double.parseDouble((String)geno.getExtendedAttribute(VcfFormatAttributes.BREAKPOINT_QUAL.attribute(), "0"));
		double bequal = Double.parseDouble((String)geno.getExtendedAttribute(VcfFormatAttributes.BREAKEND_QUAL.attribute(), "0"));
		return bpqual >= MIN_BREAKPOINT_QUAL || (bpqual == 0 & bequal >= MIN_BREAKEND_QUAL);
	}
	public static void main(String[] argv) {
        System.exit(new GeneratePonBedpe().instanceMain(argv));
    }
	public class PairBreakendSummaryWindowedSortingIterator extends WindowedSortingIterator<Pair<BreakendSummary, Integer>> {
		public PairBreakendSummaryWindowedSortingIterator(LinearGenomicCoordinate lgc, Iterator<Pair<BreakendSummary, Integer>> it, long windowSize) {
			super(it, x -> lgc.getLinearCoordinate(x.getFirst().referenceIndex, x.getFirst().start), windowSize);
		}
	}
	public static Ordering<Pair<BreakendSummary, Integer>> ByBreakendStartEnd = BreakendSummary.ByStartEnd.onResultOf(p -> p.getFirst());
}
