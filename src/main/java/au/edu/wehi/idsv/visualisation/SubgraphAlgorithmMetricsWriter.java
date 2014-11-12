package au.edu.wehi.idsv.visualisation;

import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class SubgraphAlgorithmMetricsWriter implements Closeable {
	private FileWriter fw = null;
	private BufferedWriter bw = null;
	public SubgraphAlgorithmMetricsWriter(File bedFile) throws IOException {
		try {
			fw = new FileWriter(bedFile.getAbsoluteFile());
			bw = new BufferedWriter(fw);
			//bw.append("##");
			bw.append("track gffTags=on\n");
		}
		finally {
			close();
		}
	}
	public void write(SubgraphAlgorithmMetrics metrics) throws IOException {
		bw.append(metrics.toBed());
	}
	@Override
	public void close() throws IOException {
		if (bw != null) bw.close();
		if (fw != null) fw.close();
	}
}
