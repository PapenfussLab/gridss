package au.edu.wehi.idsv.visualisation;

import java.util.List;

public interface TrackedBuffer {
	public class NamedTrackedBuffer {
		public NamedTrackedBuffer(String name, int size) {
			this.name = name;
			this.size = size;
		}
		public final String name;
		public final int size;
	}
	public void setTrackedBufferContext(String context);
	public List<NamedTrackedBuffer> currentTrackedBufferSizes();
}
