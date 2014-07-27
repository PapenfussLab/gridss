package au.edu.wehi.idsv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.SequenceUtil;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;

import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import au.edu.wehi.idsv.vcf.VcfConstants;

/**
 * Unified command line arguments for idsv
 * @author Daniel Cameron
 *
 */
public abstract class CommandLineProgram extends picard.cmdline.CommandLineProgram {
	// --- files ---
	@Option(doc="Coordinate-sorted input BAM file.",
    		shortName=StandardOptionDefinitions.INPUT_SHORT_NAME)
    public List<File> INPUT;
	@Option(doc="Coordinate-sorted tumour BAM file.",
			optional=true,
    		shortName="T")
    public List<File> INPUT_TUMOUR = Lists.newArrayList();
	@Option(doc="VCF structural variation calls.",
    		shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;
    @Option(doc="Reference used for alignment",
    		shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME)
    public File REFERENCE;
    @Option(doc="File to write recommended realignment script to",
    		optional=true)
    public File SCRIPT;
    // --- intermediate file parameters ---
    @Option(doc = "Save intermediate results into separate files for each chromosome.",
            optional = true,
            shortName = "BYCHR")
    public boolean PER_CHR = true;
    @Option(doc = "Directory to place intermediate results directories. Default location is the same directory"
    		+ " as the associated input or output file. Increases the number of intermediate files"
    		+ " but allows a greater level of parallelisation.",
            optional = true)
    public File WORKING_DIR = null;
    // --- evidence filtering parameters ---
    @Option(doc = "Minimum length of a soft-clip to be considered for analysis." +
			" Local aligners tend to produce many reads with very short soft clips.",
		optional=true)
    public int SOFT_CLIP_MIN_LENGTH = 4;
    // --- soft clip realignment parameters ---
    @Option(doc = "Minimum read MAPQ for soft clip realignment to be performed",
    		optional=true)
    public int SOFT_CLIP_REALIGNMENT_MIN_MAPQ = 5;
    @Option(doc = "Minimum sequence identity to reference. Percentage value taking a value in the range 0-100.",
    		optional=true)
    public float SOFT_CLIP_REALIGNMENT_MIN_ANCHOR_PERCENT_IDENTITY = 95;
    @Option(doc = "Local assembly algorithm used to construct breakend contigs.",
    		optional=true)
    // --- Assembly parameters ---
    public AssemblyMethod ASSEMBLY_METHOD = AssemblyMethod.DEBRUIJN_SUBGRAPH;
    @Option(doc = "k-mer used for de bruijn graph construction",
    		optional=true,
    		shortName="K")
    // --- De Bruijn assembly parameters ---
    public int ASSEMBLY_DEBRUIJN_KMER = 25;
    @Option(doc = "Maximum of base mismatches for de bruijn kmer paths to be merged",
    		optional=true)
    public int ASSEMBLY_DEBRUIJN_MAX_PATH_COLLAPSE_BASE_MISMATCHES = 2;
    @Option(doc = "Only consider bubbles for path collapse."
    		+ "Bubbles are kmer paths with a single entry and exit kmer choice",
    		optional=true)
    public boolean ASSEMBLY_DEBRUIJN_COLLAPSE_BUBBLES_ONLY = false;
    @Option(doc = "Allow reuse of reference kmers when assembling subsequent"
    		+ " contigs in an assembly iteration",
    		optional=true)
    public boolean ASSEMBLY_DEBRUIJN_ALLOW_REFERENCE_KMER_RESUSE = true;
    @Option(doc = "Maximum number of contigs per assembly iteration",
    		optional=true)
    public int ASSEMBLY_DEBRUIJN_MAX_CONTIGS_PER_ITERATION = 8;
    @Option(doc = "Number of branches consider at each kmer branch",
    		optional=true)
    public int ASSEMBLY_DEBRUIJN_SUBGRAPH_BRANCHING_FACTOR = 16;
    @Option(doc = "Number of bases (in multiples of maximum fragment size)"
    		+ "of no contributing evidence before subgraph assembly.",
    		optional=true)
    public float ASSEMBLY_DEBRUIJN_SUBGRAPH_ASSEMBLY_FRAGMENT_DELAY = 3;
    // --- Realignment parameters ---
    @Option(doc = "Minimum length of a breakend to be considered for realignment",
    		optional=true)
    public int REALIGNMENT_MIN_BREAKEND_LENGTH = 25;
    @Option(doc = "Minimum average base quality score for realignment of read soft clip",
    		optional=true)
    public float REALIGNMENT_MIN_BASE_QUALITY = 5;
    // --- output format parameters ---
	@Option(doc = "Breakends are written to VCF files as VCF v4.1 compatible breakpoints to a placeholder contig " + VcfConstants.VCF41BREAKEND_REPLACEMENT,
            optional = true,
            shortName = "VCF41")
	public boolean VCF41_COMPATIBLE = true;
	public File ensureFileExists(final File file) {
    	if (!file.exists()) {
    		throw new RuntimeException(String.format("Required input file %s is missing. Have the earlier pipeline commands executed successfully?", file));
    	}
    	return file;
    }
	public void ensureArgs() {
		IOUtil.assertFileIsReadable(REFERENCE);
		for (File f : INPUT) {
			IOUtil.assertFileIsReadable(f);
		}
		if (INPUT_TUMOUR != null) {
			for (File f : INPUT_TUMOUR) {
				IOUtil.assertFileIsReadable(f);
			}
		}
		IOUtil.assertFileIsWritable(OUTPUT);
		if (WORKING_DIR != null) {
			IOUtil.assertDirectoryIsWritable(WORKING_DIR);
		}
	}
	public void ensureDictionariesMatch() throws IOException {
		try (ReferenceSequenceFile ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE)) {
			SAMSequenceDictionary dictionary = ref.getSequenceDictionary();
			final SamReaderFactory samFactory = SamReaderFactory.makeDefault();
			for (File f : INPUT) {
				try (SamReader reader = samFactory.open(f)) {
					SequenceUtil.assertSequenceDictionariesEqual(reader.getFileHeader().getSequenceDictionary(), dictionary, f, REFERENCE);
				}
			}
			if (INPUT_TUMOUR != null) {
				for (File f : INPUT_TUMOUR) {
					try (SamReader reader = samFactory.open(f)) {
						SequenceUtil.assertSequenceDictionariesEqual(reader.getFileHeader().getSequenceDictionary(), dictionary, f, REFERENCE);
					}
				}
			}
		}
	}
	private AssemblyParameters getAssemblyParameters() {
		AssemblyParameters ap = new AssemblyParameters();
		ap.method = ASSEMBLY_METHOD;
		ap.k = ASSEMBLY_DEBRUIJN_KMER;
		ap.maxBaseMismatchForCollapse = ASSEMBLY_DEBRUIJN_MAX_PATH_COLLAPSE_BASE_MISMATCHES;
		ap.collapseBubblesOnly = ASSEMBLY_DEBRUIJN_COLLAPSE_BUBBLES_ONLY;
		ap.allReferenceKmerReuse = ASSEMBLY_DEBRUIJN_ALLOW_REFERENCE_KMER_RESUSE;
		ap.maxContigsPerAssembly = ASSEMBLY_DEBRUIJN_MAX_CONTIGS_PER_ITERATION;
		ap.subgraphAssemblyTraversalMaximumBranchingFactor = ASSEMBLY_DEBRUIJN_SUBGRAPH_BRANCHING_FACTOR;
		ap.subgraphAssemblyMargin = ASSEMBLY_DEBRUIJN_SUBGRAPH_ASSEMBLY_FRAGMENT_DELAY;
		return ap;
	}
	private SoftClipParameters getSoftClipParameters() {
		SoftClipParameters scp = new SoftClipParameters();
	    scp.minLength = SOFT_CLIP_MIN_LENGTH;
	    scp.minReadMapq = SOFT_CLIP_REALIGNMENT_MIN_MAPQ;
	    scp.minAnchorIdentity = SOFT_CLIP_REALIGNMENT_MIN_ANCHOR_PERCENT_IDENTITY;
	    return scp;
	}
	private RealignmentParameters getRealignmentParameters() {
		RealignmentParameters rp = new RealignmentParameters();
		rp.minLength = REALIGNMENT_MIN_BREAKEND_LENGTH;
		rp.minAverageQual = REALIGNMENT_MIN_BASE_QUALITY;
		return rp;
	}
	private FileSystemContext getFileSystemContext() {
		return new FileSystemContext(TMP_DIR.get(0), WORKING_DIR, MAX_RECORDS_IN_RAM);
	}
	public void close() throws IOException {
		if (processContext != null) {
			processContext.close();
			processContext = null;
		}
	}
	private ProcessingContext processContext = null;
	public ProcessingContext getContext() throws IOException {
		if (processContext == null) {
	    	processContext = new ProcessingContext(
	    			getFileSystemContext(),
	    			getDefaultHeaders(),
	    			getSoftClipParameters(),
	    			getAssemblyParameters(),
	    			getRealignmentParameters(),
	    			REFERENCE, PER_CHR, VCF41_COMPATIBLE);
		}
    	return processContext;
	}
	public static String getRealignmentScript(Iterable<EvidenceSource> it) {
    	StringBuilder sb = new StringBuilder();
    	for (EvidenceSource source : it) {
    		if (!source.isRealignmentComplete()) {
    			sb.append(source.getRealignmentScript());
    		}
    	}
    	return sb.toString();
    }
}
