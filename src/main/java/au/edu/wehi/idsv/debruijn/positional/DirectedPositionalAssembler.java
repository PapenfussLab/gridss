package au.edu.wehi.idsv.debruijn.positional;

import java.io.File;
import java.util.Iterator;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

import com.google.common.collect.PeekingIterator;

import au.edu.wehi.idsv.AssemblyEvidenceSource;
import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.DirectedEvidence;
import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMRecordAssemblyEvidence;
import au.edu.wehi.idsv.util.DuplicatingIterable;
import htsjdk.samtools.util.Log;

/**
 * Performs independent positional assembly for forward
 * and backward assembly.
 * @author Daniel Cameron
 *
 */
public class DirectedPositionalAssembler implements Iterator<SAMRecordAssemblyEvidence> {
	private static final Log log = Log.getInstance(DirectedPositionalAssembler.class);
	private static final int EVIDENCE_BUFFER_SIZE = 256;
	private static final int OUTPUT_BUFFER_SIZE = 16;
	private static final Object END_OF_STREAM = new Object();
	private final DuplicatingIterable<DirectedEvidence> forkedIterator;
	/**
	 * Assembly output buffer.
	 */
	private final BlockingQueue<Object> buffer = new ArrayBlockingQueue<Object>(OUTPUT_BUFFER_SIZE);
	private final ProcessingContext context;
	private final AssemblyEvidenceSource source;
	private final AssemblyThread forward;
	private final AssemblyThread backward;
	private volatile Exception ex;
	private volatile File forwardExport = null;
	private volatile File backwardExport = null;
	/**
	 * Count of threads that have finished calling assemblies
	 */
	private int completedCount = 0;
	private Object outputAssembly = null;
	public DirectedPositionalAssembler(ProcessingContext context, AssemblyEvidenceSource source, Iterator<DirectedEvidence> it) {
		this.context = context;
		this.source = source;
		this.forkedIterator = new DuplicatingIterable<DirectedEvidence>(2, it, EVIDENCE_BUFFER_SIZE);
		forward = new AssemblyThread(forkedIterator.iterator(), BreakendDirection.Forward);
		backward = new AssemblyThread(forkedIterator.iterator(), BreakendDirection.Backward);
		forward.start();
		backward.start();
	}
	@Override
	public boolean hasNext() {
		ensureOutput();
		return outputAssembly != null;
	}
	@Override
	public SAMRecordAssemblyEvidence next() {
		ensureOutput();
		SAMRecordAssemblyEvidence out = (SAMRecordAssemblyEvidence)outputAssembly;
		outputAssembly = null;
		return out;
	}
	private void ensureOutput() {
		if (completedCount == 2 || outputAssembly != null) return;
		do {
			try {
				outputAssembly = buffer.take();
			} catch (InterruptedException e) {
				log.error(e);
			}
			if (ex != null) {
				ex.printStackTrace();
				log.error(ex);
				throw new RuntimeException(ex);
			}
			if (outputAssembly == END_OF_STREAM) {
				completedCount++;
				outputAssembly = null;
			}
		} while (completedCount < 2 && outputAssembly == null);
	}
	public void exportGraph(File forward, File backward) {
		forwardExport = forward;
		backwardExport = backward;
	}
	private class AssemblyThread extends Thread {
		private final PositionalAssembler assembler;
		private final BreakendDirection direction;
		public AssemblyThread(final PeekingIterator<DirectedEvidence> it, final BreakendDirection direction) {
			if (it.hasNext()) {
				int referenceIndex = it.peek().getBreakendSummary().referenceIndex; 
				String chr = context.getDictionary().getSequence(referenceIndex).getSequenceName();
				setName("Assembly" + direction.toString() + chr);
			}
			this.direction = direction;
			this.assembler = new PositionalAssembler(context, source, it, direction);
		}
		@Override
		public void run() {
			try {
				while (assembler.hasNext()) {
					SAMRecordAssemblyEvidence assembly = assembler.next();
					if (forwardExport != null && direction == BreakendDirection.Forward) {
						forwardExport = null;
					} else if (backwardExport != null && direction == BreakendDirection.Backward) {
						backwardExport = null;
					}
					buffer.put(assembly);
				}
				buffer.put(END_OF_STREAM);
			} catch (Exception e) {
				ex = e;
				try {
					buffer.put(END_OF_STREAM);
				} catch (InterruptedException e1) {
					log.error(e1);
					e1.printStackTrace();
				}
			}
		}
	}
}
