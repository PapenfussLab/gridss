package au.edu.wehi.idsv;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;

/**
 * Writes the evidence stream to a file
 * @author Daniel Cameron
 *
 */
public class EvidenceToCsv {
	private final PrintStream stream;
	public EvidenceToCsv(File file) {
		try {
			this.stream = new PrintStream(new BufferedOutputStream(new FileOutputStream(file)));
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e);
		}
		writeHeader();
	}
	private void writeHeader() {
		writeCallContextHeader();
		writeDirectedEvidenceHeader();
		writeDirectedBreakpointHeader();
		stream.println();
	}
	private void writeCallContextHeader() {
		for (String s : new String[] { "callID", "assembly", "assemblyRemapped", "assemblyElevated"}) {
			stream.print(s);
			stream.print(',');
		}
	}
	private void writeDirectedEvidenceHeader() {
		for (String s : new String[] {
				"breakendQual",
				"breakpointSummary",
				//"breakendSequence",
				//"breakendQuality",
				"evidenceID",
				"localMapq",
				"exact",
				"class",
				}) {
			stream.print(s);
			stream.print(',');
		}
	}
	private void writeDirectedBreakpointHeader() {
		for (String s : new String[] {
				"breakpointQual",
				//"breakendSummary",
				"remoteMapq",
				"untemplatedSequence",
				}) {
			stream.print(s);
			stream.print(',');
		}
	}
	public void writeEvidence(DirectedEvidence evidence, VariantContextDirectedEvidence call) {
		writeSingleEvidence(evidence, null, call);
//		if (evidence instanceof AssemblyEvidence) {
//			AssemblyEvidence ass = (AssemblyEvidence)evidence;
//			for (DirectedEvidence e : ass.getEvidence()) {
//				writeSingleEvidence(e, ass, call);
//			}
//		}
	}
	private void writeSingleEvidence(DirectedEvidence evidence, AssemblyEvidence containingAssembly, VariantContextDirectedEvidence call) {
		writeCallContext(call, containingAssembly);
		writeDirectedEvidence(evidence);
		writeDirectedBreakpoint(evidence);
		stream.println();
	}
	private void writeCallContext(VariantContextDirectedEvidence call, AssemblyEvidence containingAssembly) {
		if (call != null) stream.print(call.getID());
		stream.print(',');
		if (containingAssembly != null) stream.print(containingAssembly.getEvidenceID());
		stream.print(',');
		if (containingAssembly != null) stream.print(call instanceof DirectedBreakpoint && containingAssembly instanceof DirectedBreakpoint && !call.getBreakendSummary().overlaps(containingAssembly.getBreakendSummary()));
		stream.print(',');
		if (containingAssembly != null) stream.print(!(call instanceof DirectedBreakpoint) && containingAssembly instanceof DirectedBreakpoint);
		stream.print(',');
	}
	private void writeDirectedEvidence(DirectedEvidence e) {
		stream.print(e.getBreakendQual());
		stream.print(',');
		stream.print(e.getBreakendSummary());
		stream.print(',');
		stream.print(e.getEvidenceID());
		stream.print(',');
		stream.print(e.getLocalMapq());
		stream.print(',');
		stream.print(e.isBreakendExact());
		stream.print(',');
		stream.print(e.getClass().getSimpleName());
		stream.print(',');
	}
	private void writeDirectedBreakpoint(DirectedEvidence evidence) {
		if (evidence instanceof DirectedBreakpoint) {
			DirectedBreakpoint bp = (DirectedBreakpoint)evidence;
			stream.print(bp.getBreakpointQual());
			stream.print(',');
			stream.print(bp.getRemoteMapq());
			stream.print(',');
			stream.print(bp.getUntemplatedSequence());
			stream.print(',');
		} else {
			for (int i = 0; i < 7; i++) {
				stream.print(',');
			}
		}
	}
}
