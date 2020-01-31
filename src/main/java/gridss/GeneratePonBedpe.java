package gridss;

import au.edu.wehi.idsv.*;
import au.edu.wehi.idsv.bed.BedpeWriter;
import au.edu.wehi.idsv.configuration.GridssConfiguration;
import au.edu.wehi.idsv.util.AsyncBufferedIterator;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import au.edu.wehi.idsv.vcf.VcfFormatAttributes;
import au.edu.wehi.idsv.vcf.VcfSvConstants;
import com.google.common.collect.Iterators;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;
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
    @Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output BEDPE", optional=false)
    public File OUTPUT;
	@Argument(shortName="SBO", doc="Output BED", optional=false)
	public File SINGLE_BREAKEND_OUTPUT;
	@Argument(shortName="NO", doc="VCF index of normal", optional=true)
    public int NORMAL_ORDINAL = 0;
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
	protected boolean requiresReference() {
		return true;
	}
    @Override
	protected int doWork() {
		try {
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
			Iterator<BreakendSummary> mergedIt = filteredMerge(pc, INPUT);
			BedpeMergingCounter pe = new BedpeMergingCounter();
			BedMergingCounter se = new BedMergingCounter(true);
			BedpeWriter writer = new BedpeWriter(pc.getDictionary(), OUTPUT);
			BufferedWriter sewriter = Files.newBufferedWriter(SINGLE_BREAKEND_OUTPUT.toPath(), StandardCharsets.US_ASCII);
			while(mergedIt.hasNext()) {
				BreakendSummary bs = mergedIt.next();
				if (bs instanceof BreakpointSummary) {
					writeBedpe(pe.process((BreakpointSummary)bs), writer);
				} else {
					writeBed(pc.getReference().getSequenceDictionary(), sewriter, se.process(bs));
				}
			}
			writeBedpe(pe.finish(), writer);
			writeBed(pc.getReference().getSequenceDictionary(), sewriter, se.finish());
			writer.close();
			sewriter.close();
		} catch (IOException e) {
			log.error(e);
			return 1;
		} catch (InterruptedException e) {
		}
		log.error("Imprecise call inclusion not recommended due to overly aggressive PON matching.");
		return 0;
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

	private Iterator<BreakendSummary> filteredMerge(GenomicProcessingContext pc, List<File> file) {
    	List<CloseableIterator<BreakendSummary>> fileIt = new ArrayList<>();
    	for (File f : file) {
    		fileIt.add(getFilteredIterator(pc, f));
		}
		// allocate VCF files to worker threads
		final AtomicInteger counter = new AtomicInteger(0);
		final int size = WORKER_THREADS;
		final List<Iterator<BreakendSummary>> partitioned = fileIt.stream()
				.collect(Collectors.groupingBy(it -> counter.getAndIncrement() % size))
				.values()
				.stream()
				.map(list -> list.size() == 1 ? list.get(0) : Iterators.mergeSorted(list, BreakendSummary.ByStartEnd))
				.map(it -> new AsyncBufferedIterator<>(it,"AsyncVCF"))
				.collect(Collectors.toList());
		Iterator<BreakendSummary> mergedIt = new AsyncBufferedIterator<>(Iterators.mergeSorted(partitioned, BreakendSummary.ByStartEnd), "Merged VCF reader");
		return mergedIt;
	}
	private CloseableIterator<BreakendSummary> getFilteredIterator(GenomicProcessingContext pc, File file) {
		VCFFileReader vcfReader = new VCFFileReader(file, false);
		CloseableIterator<VariantContext> it = vcfReader.iterator();
		Iterator<BreakendSummary> idsvIt = Iterators.transform(it, variant -> getBreakendSummary(pc, variant));
		Iterator<BreakendSummary> nonnullIt = Iterators.filter(idsvIt, variant -> variant != null);
		Iterator<BreakendSummary> bpit = new BreakendSummaryWindowedSortingIterator<>(pc, MAX_BREAKPOINT_HOMOLOGY_LENGTH, nonnullIt);
		return new AutoClosingIterator<>(bpit, vcfReader, it);
	}
	private BreakendSummary getBreakendSummary(GenomicProcessingContext pc, VariantContext variant) {
		IdsvVariantContext vc = IdsvVariantContext.create(pc, null, variant);
		if (!INCLUDE_IMPRECISE_CALLS && vc.hasAttribute(VcfSvConstants.IMPRECISE_KEY)) {
			return null;
		}
		if (!passesNormalFilter(vc)) {
			return null;
		}
		if (!(vc instanceof DirectedEvidence)) {
			return null;
		}
		return ((DirectedEvidence)vc).getBreakendSummary();
	}

	private boolean passesNormalFilter(IdsvVariantContext vc) {
		Genotype geno = vc.getGenotype(NORMAL_ORDINAL);
		if (geno == null) {
			throw new RuntimeException(String.format("Missing normal ordinal for ", vc.getID()));
		}
		double bpqual = Double.parseDouble((String)geno.getExtendedAttribute(VcfFormatAttributes.BREAKPOINT_QUAL.attribute(), "0"));
		double bequal = Double.parseDouble((String)geno.getExtendedAttribute(VcfFormatAttributes.BREAKEND_QUAL.attribute(), "0"));
		return bpqual >= MIN_BREAKPOINT_QUAL || (bpqual == 0 & bequal >= MIN_BREAKEND_QUAL);
	}
	public static void main(String[] argv) {
        System.exit(new GeneratePonBedpe().instanceMain(argv));
    }
}
