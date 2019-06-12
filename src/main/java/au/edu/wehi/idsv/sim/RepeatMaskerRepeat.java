package au.edu.wehi.idsv.sim;

import org.apache.commons.lang3.StringUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Represents a RepeatMasker repeat
 * 
 * TODO: proper parsing (eg: bean with flatworm)
 * 
 * @author Daniel Cameron
 *
 */
public class RepeatMaskerRepeat {
	public String chr;
	public int begin;
	public int end;
	public String repeat;
	public String classfamily;
	public RepeatMaskerRepeat(String line) {
		String[] fields = StringUtils.split(line);
		// 0      = Smith-Waterman score of the match, usually complexity adjusted
		//        The SW scores are not always directly comparable. Sometimes
		//        the complexity adjustment has been turned off, and a variety of
		//        scoring-matrices are used.
		// 1      % substitutions in matching region compared to the consensus
		// 2	  % of bases opposite a gap in the query sequence (deleted bp)
		// 3     = % of bases opposite a gap in the repeat consensus (inserted bp)
		// 4 = name of query sequence
		this.chr = fields[4].intern();
		// 5 = starting position of match in query sequence
		this.begin = Integer.parseInt(fields[5]);
		// 6 = ending position of match in query sequence
		this.end = Integer.parseInt(fields[6]);
		// 7 = no. of bases in query sequence past the ending position of match
		// 8 = match is with the Complement of the consensus sequence in the database
		// 9 = name of the matching interspersed repeat
		repeat = fields[9].intern();
		// 10 = the class of the repeat, in this case a DNA transposon 
		//		            fossil of the MER2 group (see below for list and references)
		classfamily = fields[10].intern();
		// 11 = no. of bases in (complement of) the repeat consensus sequence 
		//		            prior to beginning of the match (so 0 means that the match extended 
		//		            all the way to the end of the repeat consensus sequence)
		// 12 = starting position of match in database sequence (using top-strand numbering)
		// 13 = ending position of match in database sequence
	}
	public static List<RepeatMaskerRepeat> loadRepeatMaskerOutput(File rpOutput, String chr, String classFamily) throws IOException {
		List<RepeatMaskerRepeat> list = new ArrayList<RepeatMaskerRepeat>();
		BufferedReader br = new BufferedReader(new FileReader(rpOutput));
		// skip header lines
		br.readLine();
		br.readLine();
		br.readLine();
		for (String line = br.readLine(); line != null; line = br.readLine()) {
			RepeatMaskerRepeat r = new RepeatMaskerRepeat(line);
			if (classFamily == null || classFamily.equals(r.classfamily)) {
				if (chr == null || chr.equals(r.chr)) {
					list.add(r);
				}
			}
		}
		br.close();
		return list;
	}
}
