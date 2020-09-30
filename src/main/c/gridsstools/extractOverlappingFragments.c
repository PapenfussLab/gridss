#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <zlib.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>
#include <htslib/thread_pool.h>
#include <htslib/kseq.h>

KSTREAM_INIT(gzFile, gzread, 16384)

#define INITIAL_SIZE 1024

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   gridsstools extractOverlappingFragments [options] <region.bed> <input.bam>\n\n");
	fprintf(stderr, "Options: -o FILE  BAM output\n");
	fprintf(stderr, "         -m INT   Padding margin (Default: 2000)\n");
	fprintf(stderr, "         -@ INT   Number of threads to use\n");
	fprintf(stderr, "\n");
	return EXIT_FAILURE;
}
typedef struct interval_array {
	hts_pair_pos_t *intervals;
	int size;
	int capacity;
} interval_array;

int comp_interval_by_start(const void* a, const void* b) {
	if (((hts_pair_pos_t*)a)->beg < ((hts_pair_pos_t*)b)->beg) return -1;
	if (((hts_pair_pos_t*)a)->beg == ((hts_pair_pos_t*)b)->beg) return 0;
	return 1;
}
int comp_interval_overlaps(const void* a, const void* b) {
	hts_pair_pos_t *left = (hts_pair_pos_t*)a;
	hts_pair_pos_t *right = (hts_pair_pos_t*)b;
	if (left->end < right->beg) return -1;
	if (left->beg > right->end) return 1;
	return 0;
}
uint32_t* parse_cigar(char *str, size_t* n_cigar) {
	int n_ops = 0;
	for (char* p = str; isdigit((unsigned char)*p) || bam_cigar_table[(unsigned char)*p] >= 0; p++) {
		if (!isdigit((unsigned char)*p)) {
			n_ops++;
		}
	}
	if (n_ops == 0) return NULL;
	uint32_t *cigar = malloc(n_ops * sizeof(uint32_t));
	for (int i = 0; i < n_ops; i++) {
		int len = 0;
		while (isdigit((unsigned char)*str)) {
			len *= 10;
			len += *str - '0';
			str++;
		}
		int op = bam_cigar_table[(unsigned char)*str++];
		if (op < 0) return NULL;
		cigar[i] = (len << BAM_CIGAR_SHIFT) | op;
	}
	*n_cigar = n_ops;
	return cigar;
}
// Sorted bed intervals indexed by reference_index
interval_array* bed_to_lookup(char *filename, bam_hdr_t *header) {
	kstring_t str = { 0, 0, NULL };
	kstream_t *ks = NULL;
	gzFile fp_bed = gzopen(filename, "rb");
	if (fp_bed == NULL) {
		fprintf(stderr, "Unable to open input file %s.\n", filename);
		goto bed_to_lookup_error;
	}
	int nref = sam_hdr_nref(header);
	// initialise lookup
	interval_array *lookup = malloc(sizeof(interval_array) * nref);
	for (int i = 0; i < nref; i++) {
		lookup[i].intervals = malloc(sizeof(hts_pair_pos_t) * INITIAL_SIZE);
		lookup[i].size = 0;
		lookup[i].capacity = INITIAL_SIZE;
	}
	// read BED file
	int dret;
	int lineno = 0;
	ks = ks_init(fp_bed);
	while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
		lineno++;
		char *p, *q;
		int tid, pos, num = 0;
		int64_t beg = 0, end = 0;
		if (str.l == 0 || *str.s == '#') continue; /* empty or comment line */
		if (strncmp(str.s, "track ", 6) == 0) continue;
		if (strncmp(str.s, "browser ", 8) == 0) continue;
		for (p = q = str.s; *p && !isspace(*p); ++p);
		if (*p == 0) {
			fprintf(stderr, "Error parsing line %d of bed file.\n", lineno);
			goto bed_to_lookup_error;
		}
		char c = *p;
		*p = 0; tid = sam_hdr_name2tid(header, q); *p = c;
		if (tid < 0) {
			fprintf(stderr, "Error parsing line %d of bed file: chromosome not found in bam file.\n", lineno);
			goto bed_to_lookup_error;
		}
		num = sscanf(p + 1, "%"SCNd64" %"SCNd64, &beg, &end);
		if (num < 2 || end < beg) {
			fprintf(stderr, "Error parsing line %d of bed file.\n", lineno);
			goto bed_to_lookup_error;
		}
		if (lookup[tid].size == lookup[tid].capacity) {
			lookup[tid].intervals = realloc(lookup[tid].intervals, sizeof(hts_pair_pos_t) * lookup[tid].capacity * 2);
			lookup[tid].capacity *= 2;
		}
		lookup[tid].intervals[lookup[tid].size].beg = beg;
		lookup[tid].intervals[lookup[tid].size].end = end;
		lookup[tid].size++;
	}
	bed_to_lookup_cleanup:
	if (ks) ks_destroy(ks);
	if (fp_bed) gzclose(fp_bed);
	if (lookup) {
		for (int i = 0; i < nref; i++) {
			lookup[i].intervals = realloc(lookup[i].intervals, sizeof(hts_pair_pos_t) * lookup[i].size + 1);
			lookup[i].capacity = lookup[i].size + 1;
			qsort(lookup[i].intervals, lookup[i].size, sizeof(hts_pair_pos_t), comp_interval_by_start);
		}
	}
	return lookup;
	bed_to_lookup_error:
	// TODO cleanup
	for (int i = 0; i < nref; i++) {
		free(lookup[i].intervals);
	}
	free(lookup);
	lookup = NULL;
	goto bed_to_lookup_cleanup;
}
int overlaps(interval_array *lookup, int32_t tid, hts_pos_t start, hts_pos_t end) {
	if (tid < 0) return 0;
	hts_pair_pos_t key = { start, end };
	return bsearch(&key, lookup[tid].intervals, sizeof(hts_pair_pos_t), lookup[tid].size, comp_interval_overlaps) != NULL;
}
int any_overlap(interval_array *lookup, bam_hdr_t *header, bam1_t *record) {
	if (!(record->core.flag & BAM_FUNMAP)) {
		if (overlaps(lookup, record->core.tid, record->core.pos, bam_endpos(record) - 1)) return 1;
	}
	if (!(record->core.flag & BAM_FMUNMAP)) {
		hts_pos_t start = record->core.mpos;
		hts_pos_t end = start;
		// check for mate cigar
		kstring_t mc = KS_INITIALIZE;
		if (bam_aux_get_str(record, "MC", &mc) == 1) {
			size_t n_cigar;
			uint32_t* cigar = parse_cigar(mc.s, &n_cigar);
			if (cigar != NULL) {
				end += bam_cigar2rlen(n_cigar, cigar);
			}
			ks_free(&mc);
		}
		if (end == start) {
			// assume mate is same length as this read
			end = start + record->core.l_qseq - 1;
		}
		if (overlaps(lookup, record->core.mtid, start, end)) return 1;
	}
	kstring_t sa = KS_INITIALIZE;
	if (bam_aux_get_str(record, "SA", &sa) == 1) {
		// Iterate over SA tag alignments
		char* s = sa.s;
		while (*s != '\0') {
			// iterate through SA tag
			// we temporarily swap out the trailing seperator to make string
			// parsing of each field easier
			char *fieldstr;
			char end_of_field;
			// rname
			fieldstr = s;
			while (*s != '\0' && *s != ',' && *s != ';') s++;
			end_of_field = *s;
			*s = '\0';
			int tid = sam_hdr_name2tid(header, fieldstr);
			*s = end_of_field;
			if (*s != '\0') s++;
			// pos
			fieldstr = s;
			while (*s != '\0' && *s != ',' && *s != ';') s++;
			end_of_field = *s;
			*s = '\0';
			hts_pos_t pos = atol(fieldstr);
			*s = end_of_field;
			if (*s != '\0') s++;
			// strand
			while (*s != '\0' && *s != ',' && *s != ';') s++;
			if (*s != '\0') s++;
			// CIGAR
			fieldstr = s;
			while (*s != '\0' && *s != ',' && *s != ';') s++;
			end_of_field = *s;
			*s = '\0';
			int rlen = 0;
			size_t n_cigar = 0;
			uint32_t* cigar = parse_cigar(fieldstr, &n_cigar);
			if (cigar != NULL) {
				rlen = bam_cigar2rlen(n_cigar, cigar);
				free(cigar);
			}
			*s = end_of_field;
			if (*s != '\0') s++;
			// mapQ
			while (*s != '\0' && *s != ',' && *s != ';') s++;
			if (*s != '\0') s++;
			// NM
			while (*s != '\0' && *s != ',' && *s != ';') s++;
			if (*s != '\0') s++;
			// check for overlap
			if (tid > 0 && pos > 0 && rlen > 0) {
				if (overlaps(lookup, tid, pos, pos + rlen)) return 1;
			}
		}
		ks_free(&sa);
	}
	return 0;
}

int main_extractOverlappingFragments(int argc, char *argv[]) {
	char *out_filename = NULL;
	int threads = 1;
	int margin = 2000;
	// arg parsing
	int c;
	char *strargs = stringify_argv(argc, argv);
	while ((c = getopt(argc, argv, "o:@:m:")) >= 0) {
		switch (c) {
			case 'o': out_filename = strdup(optarg); break;
			case '@': threads = strtol(optarg, 0, 0); break;
			case 'm': margin = strtol(optarg, 0, 0); break;
			default: return usage();
		}
	}
	if (argc != optind + 2) {
		return usage();
	}
	htsFile *fp_in = NULL;
	htsFile *fp_out = NULL;
	bam_hdr_t *hdr = NULL;
	bam1_t *record = bam_init1();
	interval_array *lookup = NULL;
	int status = EXIT_SUCCESS;
	if (!(fp_in = hts_open(argv[optind + 1], "r"))) {
		fprintf(stderr, "Unable to open input file %s.\n", argv[optind + 1]);
		goto error;
	}
	htsThreadPool p = {NULL, 0};
	if (threads > 0) {
		if (!(p.pool = hts_tpool_init(threads))) {
			fprintf(stderr, "Error creating thread pool\n");
			goto error;
		}
		hts_set_opt(fp_out, HTS_OPT_THREAD_POOL, &p);
	}
	if (!(hdr = sam_hdr_read(fp_in))) {
		fprintf(stderr, "Unable to read SAM header.\n");
		goto error;
	}
	sam_hdr_add_pg(hdr, "gridsstools", "CL", strargs);
	if (!(lookup = bed_to_lookup(argv[optind], hdr))) {
		fprintf(stderr, "Unable parse target region BED.\n");
		goto error;
	}
	if (!(fp_out = hts_open(out_filename, "w"))) {
		fprintf(stderr, "Unable to open output file %s.\n", out_filename);
		goto error;
	}
	if (p.pool) {
		hts_set_opt(fp_out, HTS_OPT_THREAD_POOL, &p);
	}
	if (sam_hdr_write(fp_out, hdr) != 0) {
		fprintf(stderr, "Unable to write SAM header.\n");
		goto error;
	}
	int r;
	while ((r = sam_read1(fp_in, hdr, record)) >= 0) {
		if (any_overlap(lookup, hdr, record)) {
			if (sam_write1(fp_out, hdr, record) < 0) {
				fprintf(stderr, "Unable to write to output file %s.\n", out_filename);
				goto error;
			}
		}
	}
	if (r < -1) {
		fprintf(stderr, "Truncated file.\n");
		goto error;
	}
cleanup:
	if (lookup) {
		for (int i = 0; i < sam_hdr_nref(hdr); i++) {
			if (lookup[i].intervals) {
				free(lookup[i].intervals);
			}
		}
		free(lookup);
	}
	if (strargs) free(strargs);
	if (fp_out) hts_close(fp_out);
	if (fp_in) hts_close(fp_in);
	if (record) bam_destroy1(record);
	if (hdr) bam_hdr_destroy(hdr);
	return status;
error:
	status = EXIT_FAILURE;
	perror("extractOverlappingFragments failed");
	goto cleanup;
}