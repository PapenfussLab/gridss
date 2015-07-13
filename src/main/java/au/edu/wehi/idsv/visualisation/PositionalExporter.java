package au.edu.wehi.idsv.visualisation;

import htsjdk.samtools.util.CloserUtil;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.List;

import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.debruijn.positional.KmerPathNode;

/**
 * Exports a positional graph to a collection of FASTG files
 * @author cameron.d
 *
 */
public class PositionalExporter {
	/**
	 * Export loaded graph
	 * @throws IOException 
	 */
	public static void exportFastg(File file, int k, Collection<KmerPathNode> nodes) throws IOException {
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(file));
			writer.append("#FASTG:begin:version=1.0;\n");
			for (KmerPathNode n : nodes) {
				writer.append('>');
				writer.append(id(n, k));
				List<KmerPathNode> nextList = n.next();
				if (nextList.size() > 0) {
					writer.append(':');
					writer.append(id(nextList.get(0), k));
					for (int i = 1; i < nextList.size(); i++) {
						writer.append(',');
						writer.append(id(nextList.get(i), k));
					}
				}
				writer.append(":start=");
				writer.append(Integer.toString(n.firstStart()));
				writer.append(",end=");
				writer.append(Integer.toString(n.firstEnd()));
				writer.append(",weight=");
				writer.append(Integer.toString(n.weight()));
				writer.append(",reference=");
				writer.append(n.isReference() ? '1' : '0');
				writer.append(";\n");
				writer.append(new String(KmerEncodingHelper.baseCalls(n.pathKmers(), k)));
				writer.append('\n');
			}
			writer.append("#FASTG:end;\n");
		} finally {
			CloserUtil.close(writer);
		}
	}
	public static String id(KmerPathNode node, int k) {
		String str = String.format("%s_%d", new String(KmerEncodingHelper.encodedToPicardBases(k, node.firstKmer())), node.firstStart());
		if (node.firstStart() < 0) {
			str = str.replace('-', '_');
		}
		return str;
	}
	/**
	 * Exports to graphviz dot format.
	 * Attempts to be consistent with ABySS attribute names but these are poorly
	 * documented in ABySS.
	 * @param file
	 * @param k
	 * @param nodes
	 * @throws IOException
	 */
	public static void exportDot(File file, int k, Collection<KmerPathNode> nodes) throws IOException {
		// d distance between contigs (negative indicates overlap)
		// l length
		// C kmer coverage
		// writer.append(String.format("edge [d=-%d]\n", k-1));
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(file));
			writer.append("digraph {\n");
			writer.append("rankdir=LR;");
			writer.append(String.format("graph [k=%d]\n", k));
			for (KmerPathNode n : nodes) {
				writer.append(id(n, k));
				writer.append(String.format(" [s=%d,e=%d,wid=%d,w=%d,l=%d,r=%s,alt=%d,s=\"%s\"];\n",
						n.firstStart(), n.firstEnd(), n.width(), n.weight(), n.length(), n.isReference() ? "true" : "false", n.collapsedKmers().size(), new String(KmerEncodingHelper.baseCalls(n.pathKmers(), k))));
				for (KmerPathNode next : n.next()) {
					writer.append(id(n,k));
					writer.append(" -> ");
					writer.append(id(next,k));
					writer.append(";");
				}
				writer.append(new String(KmerEncodingHelper.baseCalls(n.pathKmers(), k)));
				
			}
			writer.append("}\n");
		} finally {
			CloserUtil.close(writer);
		}
	}
	public static void exportfasta(File file, int k, Collection<KmerPathNode> nodes) throws IOException {
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(file));
			for (KmerPathNode n : nodes) {
				writer.append('>');
				writer.append(id(n, k));
				writer.append('\n');
				writer.append(new String(KmerEncodingHelper.baseCalls(n.pathKmers(), k)));
				writer.append('\n');
			}
		} finally {
			CloserUtil.close(writer);
		}
	}
}
