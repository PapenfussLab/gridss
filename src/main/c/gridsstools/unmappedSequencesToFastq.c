#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <unistd.h>
#include <ctype.h>
#include "fastq.h"
#include "htslib/kbitset.h"

#define MAX_CONTIG_NAME_LENGTH 254

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   gridsstools unmappedSequencesToFastq [options] [input.bam] ...\n\n");
	fprintf(stderr, "Options: -o FILE  output to file\n");
	fprintf(stderr, "         -n FILE  Read names of fragments with any reference alignment.\n");
	fprintf(stderr, "         -c FILE  File containing contigs to treat as unmapped. One contig per line.\n");
	fprintf(stderr, "                  Useful when the reference genome contains decoy contigs.\n");
	fprintf(stderr, "         -m INT   Minimum length of sequence export. Minimum length of sequence export. [20]\n");
	fprintf(stderr, "         -x       exclude soft clipped bases.\n");
	fprintf(stderr, "         -u       Ensure exported names are unique by suffixing with '/1' or '/2'\n");
	fprintf(stderr, "         -@       Number of threads to use\n");
	fprintf(stderr, "\n");
	return EXIT_FAILURE;
}

static kbitset_t* load_unmapped_contig_lookup(sam_hdr_t *hdr, char* filename) {
	FILE *fp = fopen(filename, "r");
	if (!fp) return NULL;
	kbitset_t* lookup = kbs_init(hdr->n_targets);
	size_t buffersize = MAX_CONTIG_NAME_LENGTH + 1;
	char *buffer = malloc(sizeof(char) * buffersize);
	int line_length;
	while ((line_length = getline(&buffer, &buffersize, fp)) != -1) {
		if (buffer[line_length - 1] == '\n') {
			buffer[line_length - 1] = '\0';
		}
		int tid = sam_hdr_name2tid(hdr, buffer);
		if (tid >= 0) {
			kbs_insert(lookup, tid);
		}
	}
	fclose(fp);
	return lookup;
}

static int parse_cigar_op(
		char** str,
		int* op,
		int* len) {
	*len = 0;
	while (isdigit(**str)) {
		*len *= 10;
		*len += **str - '0';
		(*str)++;
	}
	if (len == 0) {
		return 0;
	}
	*op = bam_cigar_table[(unsigned char)**str];
	(*str)++;
	return *op != -1;
}
// Extracts the start and end clipping lengths from the chimeric alignment in the given SA tag
// str: [in/out] parameter of current position in the SA string.
// 		output value is the position after the SA chimeric alignment,
// 		or the position the parsing failure
// start_clip_length: [out] length of start clipping (hard and soft)
// end_clip_length: [out] length of end clipping (hard and soft)
// returns true if the parsing was successful
static int parse_SA_tag_clipping(
		char** str,
		char* reference_name,
		int* start_clip_length,
		int* end_clip_length,
		int* is_negative_strand) {
	// SA:Z:(rname ,pos ,strand ,CIGAR ,mapQ ,NM ;)+
	char* rname = *str;
	*start_clip_length = 0;
	*end_clip_length = 0;
	*is_negative_strand = 0;
	*str = strchr(*str, ','); if (!*str) return 0;
	int rname_length = *str - rname;
	rname_length = rname_length > MAX_CONTIG_NAME_LENGTH ? MAX_CONTIG_NAME_LENGTH : rname_length;
	strncpy(reference_name, rname, rname_length);
	reference_name[rname_length] = '\0';
	*str += 1;
	*str = strchr(*str, ','); if (!*str) return 0;
	*str += 1;
	*is_negative_strand = **str == '-';
	*str = strchr(*str, ','); if (!*str) return 0;
	*str += 1;
	// process CIGAR
	int oplen;
	int op;
	int inStart = 1;
	while (parse_cigar_op(str, &op, &oplen)) {
		if (op == BAM_CHARD_CLIP || op == BAM_CSOFT_CLIP) {
			if (inStart) {
				*start_clip_length += oplen;
			} else {
				*end_clip_length += oplen;
			}
		} else {
			inStart = 0;
		}
	}
	// advance to end of chimeric alignment record
	while (*str && **str != ';' && **str != '\0') {
		(*str)++;
	}
	if (*str && **str == ';') (*str)++;
	return 1;
}

static int process(
		BGZF *out,
		FILE* readname_fptr,
		kstring_t *linebuf,
		int min_sequence_length,
		int include_soft_clips,
		int unique_names,
		kbitset_t *unmapped_contig_lookup,
		sam_hdr_t *hdr,
		bam1_t *record) {
	char sa_reference_name_buffer[MAX_CONTIG_NAME_LENGTH + 1];
	int i;
	if (record->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) {
		// Only process primary alignments
		return 0;
	}
	int start = 0, end = 0;
	if (record->core.flag & BAM_FUNMAP) {
		// unmapped = full length
		start = 0;
		end = record->core.l_qseq;
	} else if (unmapped_contig_lookup && kbs_exists(unmapped_contig_lookup, record->core.tid)) {
		// Grab the full read for now.
		// Technically we should check if it's a split read alignment and exclude the
		// reference-aligning portion of the read
		start = 0;
		end = record->core.l_qseq;
		// TODO: reduce the exported sequence length by the length of any chimeric alignment
	} else if (include_soft_clips || (unmapped_contig_lookup && kbs_exists(unmapped_contig_lookup, record->core.tid))) {
		// check for start/end soft clip
		int start_clip_length = 0;
		int end_clip_length = 0;
		if (unmapped_contig_lookup && kbs_exists(unmapped_contig_lookup, record->core.tid)) {
			// treat this alignment as unmapped
			// if there's a chimeric match to the reference then we want to
			// exclude that sequence
			// we can do so using the same logic as the clipping logic by treating the entire read as clipped
			start_clip_length = record->core.l_qseq;
			end_clip_length = record->core.l_qseq;
		} else {
			uint32_t* cigar = bam_get_cigar(record);
			int n_cigar = record->core.n_cigar;
			for (i = 0; i < n_cigar; i++) {
				if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
					start_clip_length += bam_cigar_oplen(cigar[i]);
				} else if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
				} else {
					break;
				}
			}
			for (i = n_cigar - 1; i >= 0; i--) {
				if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
					end_clip_length += bam_cigar_oplen(cigar[i]);
				} else if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
				} else {
					break;
				}
			}
		}
		// SA:Z:(rname ,pos ,strand ,CIGAR ,mapQ ,NM ;)+
		uint8_t *satag = bam_aux_get(record, "SA");
		if (satag) {
			char *sa = bam_aux2Z(satag);
			if (sa) {
				int sa_start_clip_length = 0;
				int sa_end_clip_length = 0;
				int saNegativeStrand;
				while(parse_SA_tag_clipping(&sa, sa_reference_name_buffer, &sa_start_clip_length, &sa_end_clip_length, &saNegativeStrand)) {
					int sa_tid = sam_hdr_name2tid(hdr, sa_reference_name_buffer);
					if (sa_tid < 0 || kbs_exists(unmapped_contig_lookup, sa_tid)) {
						// split read to a viral reference contig or one to a contig that's not actually in the reference
						// ignore this SA alignment record
					} else {
						if ((saNegativeStrand && !bam_is_rev(record)) ||
								(!saNegativeStrand && bam_is_rev(record))) {
							// strand differs from our record - need to flip start and end
							int tmp = sa_start_clip_length;
							sa_start_clip_length = sa_end_clip_length;
							sa_end_clip_length = tmp;
						}
						start_clip_length = sa_start_clip_length < start_clip_length ? sa_start_clip_length : start_clip_length;
						end_clip_length = sa_end_clip_length < end_clip_length ? sa_end_clip_length : end_clip_length;
					}
				}
			}
		}
		if (start_clip_length >= min_sequence_length && start_clip_length > end_clip_length) {
			start = 0;
			end = start_clip_length;
		} else if (end_clip_length >= min_sequence_length) {
			end = record->core.l_qseq;
			start = record->core.l_qseq - end_clip_length;
		}
	}
	if (end - start >= min_sequence_length) {
		write_fastq(record, out, linebuf, start, end, unique_names);
		if (readname_fptr) {
			if (!(record->core.flag & BAM_FUNMAP) || (record->core.flag & BAM_FPAIRED && !(record->core.flag & BAM_FMUNMAP))) {
				fputs(bam_get_qname(record), readname_fptr);
				fputc('\n', readname_fptr);
			}
		}
		return 1;
	}
	return 0;
}

int main_unmappedSequencesToFastq(int argc, char *argv[]) {
	char *out_filename = "-";
	char *readname_filename = NULL;
	char *unmappedcontig_filename = NULL;
	int min_sequence_length = 20, include_soft_clips = 1, unique_names = 0, threads = 1;
	htsFile *fp_in = NULL;
	BGZF *fp_out = NULL;
	bam_hdr_t *hdr = NULL;
	bam1_t *record = bam_init1();
	FILE* readname_fptr = NULL;
	int status = EXIT_SUCCESS;
	kbitset_t *is_unmapped_contig = NULL;
	kstring_t linebuf = KS_INITIALIZE;
	int i;
	// arg parsing
	int c;
	while ((c = getopt(argc, argv, "o:n:c:m:xu@:")) >= 0) {
		switch (c) {
			case 'o': out_filename = strdup(optarg); break;
			case 'n': readname_filename = strdup(optarg); break;
			case 'c': unmappedcontig_filename = strdup(optarg); break;
			case 'm': min_sequence_length = strtol(optarg, 0, 0); break;
			case 'x': include_soft_clips = 0; break;
			case 'u': unique_names = 1; break;
			case '@': threads = strtol(optarg, 0, 0); break;
			default: return usage();
		}
	}
	if (argc < optind + 1) {
		return usage();
	}
	for (i = optind; i < argc; i++) {
		char* in_filename = argv[i];
		// open files
		if (!(fp_in = hts_open(in_filename, "r"))) {
			fprintf(stderr, "Unable to open input file %s.\n", in_filename);
			goto error;
		}
		if (hts_set_threads(fp_in, threads)) {
			fprintf(stderr, "Unable to create thread pool with %d threads\n", threads);
			goto error;
		}
		if (!(hdr = sam_hdr_read(fp_in))) {
			fprintf(stderr, "Unable to read SAM header.\n");
			goto error;
		}
		if (is_unmapped_contig) {
			// each input file could have different reference genomes
			kbs_destroy(is_unmapped_contig);
			is_unmapped_contig = NULL;
		}
		is_unmapped_contig = load_unmapped_contig_lookup(hdr, unmappedcontig_filename);
		if (!(fp_out = bgzf_open(out_filename, "wu"))) {
			fprintf(stderr, "Unable to open output file\n");
			goto error;
		}
		if (readname_filename && !(readname_fptr = fopen(readname_filename, "wb"))) {
			fprintf(stderr, "Unable to open read name output file\n");
			goto error;
		}
		// process records
		int r;
		while ((r = sam_read1(fp_in, hdr, record)) >= 0) {
			process(fp_out, readname_fptr, &linebuf, min_sequence_length, include_soft_clips, unique_names, is_unmapped_contig, hdr, record);
		}
		if (r < -1) {
			fprintf(stderr, "Truncated file.\n");
			goto error;
		}
		if (hdr) {
			bam_hdr_destroy(hdr);
			hdr = NULL;
		}
		if (fp_in) {
			hts_close(fp_in);
			fp_in = NULL;
		}
	}
cleanup:
	ks_free(&linebuf);
	if (record) bam_destroy1(record);
	if (fp_out) bgzf_close(fp_out);
	if (hdr) bam_hdr_destroy(hdr);
	if (fp_in) hts_close(fp_in);
	if (is_unmapped_contig) kbs_destroy(is_unmapped_contig);
	return status;
error:
	status = EXIT_FAILURE;
	perror("main_unmappedSequencesToFastq failed");
	goto cleanup;
}