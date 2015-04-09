package au.edu.wehi.idsv;

import java.util.concurrent.atomic.AtomicInteger;

public class SequentialIdGenerator implements VariantIdGenerator, AssemblyIdGenerator {
	private final AtomicInteger id = new AtomicInteger(0);
	private final String prefix;
	private final String suffix;
	public SequentialIdGenerator(String prefix) {
		this(prefix, "");
	}
	public SequentialIdGenerator(String prefix, String suffix) {
		this.prefix = prefix;
		this.suffix = suffix;
	}
	public String generate() {
		return String.format("%s%d%s", prefix, id.incrementAndGet(), suffix);
	}
	public String generate(BreakpointSummary breakpoint) {
		return generate();
	}
	@Override
	public String generate(BreakendSummary breakpoint, byte[] baseCalls, int startAnchoredBaseCount, int endAnchoredBaseCount) {
		return generate();
	}
}
