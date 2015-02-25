package au.edu.wehi.idsv;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import au.edu.wehi.idsv.sam.SAMRecordUtil;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

/**
 * Assembly spanning a small indel
 */
public class SmallIndelSAMRecordAssemblyEvidence extends RealignedSAMRecordAssemblyEvidence {
	
	private SAMRecord assembly;
	public SmallIndelSAMRecordAssemblyEvidence(AssemblyEvidenceSource source, SAMRecord assembly) {
		this(source, assembly, new IndelInfo(assembly));
	}
	private SmallIndelSAMRecordAssemblyEvidence(AssemblyEvidenceSource source, SAMRecord assembly, IndelInfo info) {
		super(source, info.createAnchored(), info.createRealigned());
		this.assembly = assembly;
	}
	private static class IndelInfo {
		private SAMRecord read;
		private BreakendDirection dir;
		private CigarElement op;
		private List<CigarElement> anchorCigar;
		private List<CigarElement> realignCigar;
		public IndelInfo(SAMRecord read) {
			this.read = read;
			this.dir = getBreakendDirection(read);
			calcCigar();
		}
		private void calcCigar() {
			int offset = -1;
			int opSize = 0;
			List<CigarElement> list = read.getCigar().getCigarElements();
			for (int i = 0; i < list.size(); i++) {
				CigarElement e = list.get(i);
				// Just grab the biggest one
				if (e.getLength() > opSize && (e.getOperator() == CigarOperator.DELETION || e.getOperator() == CigarOperator.INSERTION)) {
					offset = i;
					opSize = e.getLength();
					op = e;
				}
			}
			if (offset < 0) {
				throw new IllegalArgumentException(String.format("Indel not found in %s", read.getReadName()));
			}
			List<CigarElement> pre = list.subList(0, offset);
			List<CigarElement> post = list.subList(offset + 1, list.size());
			// assign to anchor/realign based on dir
			if (dir == BreakendDirection.Forward) {
				anchorCigar = pre;
				realignCigar = post;
			} else {
				anchorCigar = post;
				realignCigar = pre;
			}
		}
		/**
		 * Emulates the anchored portion as if the assembly was a split read mapping
		 * @param indelAssembly assembly spanning indel breakpoint
		 * @param 
		 * @return anchored breakend alignment
		 */
		public SAMRecord createAnchored() {
			SAMRecord r = SAMRecordUtil.clone(read);
			List<CigarElement> cigarList = new ArrayList<CigarElement>(anchorCigar);
			cigarList.add(dir == BreakendDirection.Forward ? cigarList.size() : 0, new CigarElement(r.getReadLength() - new Cigar(anchorCigar).getReadLength(), CigarOperator.SOFT_CLIP));
			Cigar cigar = new Cigar(cigarList);
			int offset = 0;
			if (dir == BreakendDirection.Backward) {
				offset = r.getCigar().getReferenceLength() - cigar.getReferenceLength();
			} 
			r.setAlignmentStart(r.getAlignmentStart() + offset);
			r.setCigar(cigar);
			return r;
		}
		/**
		 * Emulates a split read realignment as if the assembly was a split read mapping
		 * @param indelAssembly assembly spanning indel breakpoint
		 * @return breakend split read realignment
		 */
		private SAMRecord createRealigned() {
			SAMRecord r = SAMRecordUtil.clone(read);
			List<CigarElement> cigarList = new ArrayList<CigarElement>(realignCigar);
			if (op.getOperator() == CigarOperator.I) {
				cigarList.add(dir == BreakendDirection.Forward ? 0 : cigarList.size(), new CigarElement(op.getLength(), CigarOperator.SOFT_CLIP));
			}
			Cigar cigar = new Cigar(cigarList);
			int refOffset = 0, readOffset = 0;
			if (dir == BreakendDirection.Forward) {
				refOffset = r.getCigar().getReferenceLength() - cigar.getReferenceLength();
				readOffset = new Cigar(anchorCigar).getReadLength();
			}
			r.setAlignmentStart(r.getAlignmentStart() + refOffset);
			r.setCigar(cigar);
			r.setReadBases(Arrays.copyOfRange(r.getReadBases(), readOffset, readOffset + cigar.getReadLength()));
			r.setBaseQualities(Arrays.copyOfRange(r.getBaseQualities(), readOffset, readOffset + cigar.getReadLength()));
			return r;
		}
	}
	public SAMRecord getIndelSAMRecord() {
		return assembly;
	}
}
