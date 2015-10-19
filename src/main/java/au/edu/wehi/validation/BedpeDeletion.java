package au.edu.wehi.validation;

public class BedpeDeletion {
	public static final String FS = "\t";
	public final String line;
	public final String chrom1;
	public final int start1;
	public final int end1;
	public final String chrom2;
	public final int start2;
	public final int end2;
	//public final String name;
	//public final double score;
	public final int length;
	public BedpeDeletion(String line) {
		this.line = line;
		String[] records = line.split(FS); 
		this.chrom1 = records[0];
		this.start1 = Integer.parseInt(records[1]);
		this.end1 = Integer.parseInt(records[2]);
		this.chrom2 = records[3];
		this.start2 = Integer.parseInt(records[4]);
		this.end2 = Integer.parseInt(records[5]);
		//this.name = records[6];
		//this.score = Integer.parseInt(records[7]);
		this.length = Integer.parseInt(records[10]);
	}
}
