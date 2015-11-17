package au.edu.wehi.idsv.visualisation;

import htsjdk.samtools.util.CloserUtil;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import au.edu.wehi.idsv.debruijn.KmerEncodingHelper;
import au.edu.wehi.idsv.debruijn.positional.ImmutableKmerNode;
import au.edu.wehi.idsv.debruijn.positional.KmerNode;
import au.edu.wehi.idsv.debruijn.positional.KmerNodeUtil;
import au.edu.wehi.idsv.debruijn.positional.KmerPathNode;
import au.edu.wehi.idsv.debruijn.positional.KmerPathSubnode;

/**
 * Exports a positional graph to a collection of FASTG files
 * @author Daniel Cameron
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
	public static String id(KmerNode n, int k) {
		String str = String.format("%s_%d", new String(KmerEncodingHelper.encodedToPicardBases(k, n.firstKmer())), n.firstStart());
		if (n.firstStart() < 0) {
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
	 * @param nodes graph nodes
	 * @param contig graph nodes being assembled
	 * @throws IOException
	 */
	public static void exportDot(File file, int k, Collection<KmerPathNode> nodes, Collection<KmerPathSubnode> contig) throws IOException {
		// d distance between contigs (negative indicates overlap)
		// l length
		// C kmer coverage
		// writer.append(String.format("edge [d=-%d]\n", k-1));
		Set<KmerNode> lookup = new TreeSet<KmerNode>(KmerNodeUtil.ByFirstStartKmer);
		if (contig != null) {
			lookup.addAll(contig.stream().map(sn -> sn.node()).collect(Collectors.toList()));
		}
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(file));
			writer.append("digraph G {\n");
			//writer.append("rankdir=LR;");
			writer.append(String.format("	graph [k=%d];\n", k));
			// nodes
			for (KmerPathNode n : nodes) {
				writer.append('\t');
				writer.append(id(n, k));
				writer.append(String.format(" [s=%d,e=%d,wid=%d,w=%d,l=%d,r=%s,alt=%d,seq=\"%s\",contig=%s];\n",
						n.firstStart(), n.firstEnd(), n.width(), n.weight(), n.length(), n.isReference() ? "true" : "false",
								n.collapsedKmers().size(), new String(KmerEncodingHelper.baseCalls(n.pathKmers(), k)),
								lookup.contains(n) ? "true" : "false"));
			}
			// edges
			for (KmerPathNode n : nodes) {
				for (KmerPathNode next : n.next()) {
					writer.append('\t');
					writer.append(id(n,k));
					writer.append(" -> ");
					writer.append(id(next,k));
					writer.append(String.format(" [seq=\"%s\"];\n", KmerEncodingHelper.lastBaseEncodedToPicardBase(next.firstKmer())));
					writer.append(";\n");
				}
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
	/**
	 * Exports a full-size uncompressed positional de Bruijn graph to graphviz dot format
	 */
	public static void exportNodeDot(File file, int k, Collection<KmerPathNode> graph, Collection<KmerPathSubnode> contig) throws IOException {
		List<KmerNode> nodes = ImmutableKmerNode.split(graph).collect(Collectors.toList());
		Set<KmerNode> lookup = new TreeSet<KmerNode>(KmerNodeUtil.ByFirstStartKmer);
		lookup.addAll(nodes);
		TreeSet<KmerNode> contigLookup = new TreeSet<KmerNode>(KmerNodeUtil.ByFirstStartKmer);
		if (contig != null) {
			contigLookup.addAll(ImmutableKmerNode.split(contig).collect(Collectors.toList()));
		}
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(file));
			writer.append("digraph G {\n");
			writer.append(String.format("\tgraph [k=%d];\n", k));
			// nodes
			for (KmerNode n : nodes) {
				writer.append('\t');
				writer.append(id(n, k));
				writer.append(String.format(" [s=%d,w=%d,r=%s,seq=\"%s\",contig=%s];\n",
					n.firstStart(),n.weight(), n.isReference() ? "true" : "false", new String(KmerEncodingHelper.encodedToPicardBases(k, n.firstKmer())),
					contigLookup.contains(n) ? "true" : "false"));
			}
			// edges
			for (KmerNode n : nodes) {
				for (long nextkmer : KmerEncodingHelper.nextStates(k, n.lastKmer())) {
					ImmutableKmerNode next = new ImmutableKmerNode(nextkmer, n.firstEnd() + 1, n.firstEnd() + 1, false, 0);
					if (lookup.contains(next)) {
						writer.append('\t');
						writer.append(id(n,k));
						writer.append(" -> ");
						writer.append(id(next, k));
						writer.append(String.format(" [seq=\"%s\"];\n", KmerEncodingHelper.lastBaseEncodedToPicardBase(next.firstKmer())));
						writer.append(";\n");
					}
				}
			}
			writer.append("}\n");
		} finally {
			CloserUtil.close(writer);
		}
	}
}
