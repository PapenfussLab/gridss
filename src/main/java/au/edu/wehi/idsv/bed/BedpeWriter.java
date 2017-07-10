package au.edu.wehi.idsv.bed;

import java.io.BufferedOutputStream;
import java.io.Closeable;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.BreakpointSummary;
import au.edu.wehi.idsv.VariantContextDirectedBreakpoint;
import au.edu.wehi.idsv.vcf.VcfInfoAttributes;
import au.edu.wehi.idsv.vcf.VcfSvConstants;
import htsjdk.samtools.SAMSequenceDictionary;

/**
 * Outputs variant calls to BEDPE
 * @author Daniel Cameron
 *
 */
public class BedpeWriter implements Closeable {
	private boolean includeGridssVcfFields;
	private boolean includeUntemplatedSequence;
	private SAMSequenceDictionary dict;
	private OutputStream os;
	public BedpeWriter(SAMSequenceDictionary dictionary, File file) throws IOException {
		this(dictionary, new FileOutputStream(file));
	}
	public BedpeWriter(SAMSequenceDictionary dictionary, OutputStream stream) throws IOException {
		this.dict = dictionary;
		this.os = new BufferedOutputStream(stream);
	}
	public void writeHeader(boolean includeUntemplatedSequence, boolean includeGridssVcfFields) throws IOException {
		this.includeUntemplatedSequence = includeUntemplatedSequence;
		this.includeGridssVcfFields = includeGridssVcfFields;
		String str = "#" +
				"chrom1" + "\t" +
				"start1" + "\t" +
				"end1" + "\t" +
				"chrom2" + "\t" +
				"start2" + "\t" +
				"end2" + "\t" +
				"name" + "\t" +
				"score" + "\t" +
				"strand1" + "\t" +
				"strand2";
		if (includeUntemplatedSequence) {
			str += "\t" + "untemplatedSequence";
		}
		if (includeGridssVcfFields) {
			str += "\t" + VcfSvConstants.HOMOLOGY_SEQUENCE_KEY + "\t" +
					VcfInfoAttributes.BREAKPOINT_ASSEMBLY_COUNT.attribute() + "\t" +
					VcfInfoAttributes.BREAKPOINT_ASSEMBLY_COUNT_REMOTE.attribute() + "\t" +
					VcfInfoAttributes.BREAKPOINT_SPLITREAD_COUNT.attribute() + "\t" +
					VcfInfoAttributes.BREAKPOINT_READPAIR_COUNT.attribute() + "\t" +
					VcfInfoAttributes.BREAKPOINT_ASSEMBLY_QUAL.attribute() + "\t" +
					VcfInfoAttributes.BREAKPOINT_ASSEMBLY_QUAL_REMOTE.attribute() + "\t" +
					VcfInfoAttributes.BREAKPOINT_SPLITREAD_QUAL.attribute() + "\t" +
					VcfInfoAttributes.BREAKPOINT_READPAIR_QUAL.attribute() + "\t" +
					VcfInfoAttributes.REFERENCE_READ_COUNT.attribute() + "\t" +
					VcfInfoAttributes.REFERENCE_READPAIR_COUNT.attribute();
		}
		str += "\n";
		os.write(str.getBytes(StandardCharsets.UTF_8));
	}
	public void write(VariantContextDirectedBreakpoint variant) throws IOException {
		BreakpointSummary bp = variant.getBreakendSummary();
		List<String> args = new ArrayList<>();
		if (includeUntemplatedSequence) {
			args.add(variant.getHomologySequence());
		}
		if (includeGridssVcfFields) {
			args.addAll(Arrays.asList(new String[] {
				Integer.toString(variant.getBreakpointEvidenceCountLocalAssembly()),
				Integer.toString(variant.getBreakpointEvidenceCountRemoteAssembly()),
				Integer.toString(variant.getBreakpointEvidenceCountSoftClip()),
				Integer.toString(variant.getBreakpointEvidenceCountReadPair()),
				Double.toString(variant.getBreakpointEvidenceQualLocalAssembly()),
				Double.toString(variant.getBreakpointEvidenceQualRemoteAssembly()),
				Double.toString(variant.getBreakpointEvidenceQualSoftClip()),
				Double.toString(variant.getBreakpointEvidenceQualReadPair()),
				Integer.toString(variant.getReferenceReadCount()),
				Integer.toString(variant.getReferenceReadPairCount())
			}));
		}
		write(bp,
				variant.getID(),
				Double.toString(variant.getPhredScaledQual()),
				args.toArray(new String[0]));
	}
	public void write(BreakpointSummary bp, String id, String score, String... fields) throws IOException {
		StringBuilder sb = new StringBuilder();
		sb.append(dict.getSequence(bp.referenceIndex).getSequenceName());
		sb.append('\t'); sb.append(Integer.toString(bp.start - 1));
		sb.append('\t'); sb.append(Integer.toString(bp.end));
		sb.append('\t'); sb.append(dict.getSequence(bp.referenceIndex2).getSequenceName());
		sb.append('\t'); sb.append(Integer.toString(bp.start2 - 1));
		sb.append('\t'); sb.append(Integer.toString(bp.end2));
		sb.append('\t'); sb.append(id);
		sb.append('\t'); sb.append(score);
		sb.append('\t'); sb.append(bp.direction == BreakendDirection.Forward ? '+' : '-');
		sb.append('\t'); sb.append(bp.direction2 == BreakendDirection.Forward ? '+' : '-');
		for (String s : fields) {
			sb.append('\t'); sb.append(s);
		}
		sb.append('\n');
		os.write(sb.toString().getBytes(StandardCharsets.UTF_8));
	}
	@Override
	public void close() throws IOException {
		if (os != null) os.close();
		os = null;
	}
}
