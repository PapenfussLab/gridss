#include <stdio.h>
#include <string.h>
#include "config.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "1.0"
#endif

int main_extractFragmentsToFastq(int argc, char *argv[]);
int main_unmappedSequencesToFastq(int argc, char *argv[]);

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: gridsstools (C port of some GRIDSS tools)\n");
	fprintf(stderr, "Version: %s\n\n", PACKAGE_VERSION);
	fprintf(stderr, "Usage:   gridsstools <command> [options]\n\n");
	fprintf(stderr, "Command: extractFragmentsToFastq        Extracts a subset of reads to fastq\n");
	fprintf(stderr, "         unmappedSequencesToFastq       exports unmapped reads and bases to fastq\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc < 2) return usage();
	if (strcmp(argv[1], "extractFragmentsToFastq") == 0) return main_extractFragmentsToFastq(argc-1, argv+1);
	else if (strcmp(argv[1], "unmappedSequencesToFastq") == 0) return main_unmappedSequencesToFastq(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}
	return 0;
}