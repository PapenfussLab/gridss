import argparse
import pysam

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Revert supplementary UA alignments with low mapping quality')
parser.add_argument('--before', required=True, help='The file before UA')
parser.add_argument('--after', required=True, help='The file after UA')
parser.add_argument('--output', required=True, help='The output realigned file')
parser.add_argument('--min_mapping_quality', type=int, default=60, help='low mapping quality threshold')
args = parser.parse_args()

# Open BAM files - before and after UA alignments
before = pysam.AlignmentFile(args.before, "rb")
after = pysam.AlignmentFile(args.after, "rb")

# Create a new BAM file for output
output = pysam.AlignmentFile(args.output, "wb", template=after)

reads_to_remove = set()
sup_reads = set()

# added 'fetch' to make sure we iterate over all reads - otherwise it does not iterate once more

# reads with supplementary alignment
for after_read in after.fetch():
    if after_read.is_supplementary:
        sup_reads.add(after_read.query_name)

# reads that should be reverted:
# reads with supplementary alignment (flag = 2048) and has low mapping quality
# reads with decoy alignment
# reads which their primary mapping (flag = 0) quality is low, and they have supplementary alignment read
for after_read in after.fetch():
        if (after_read.mapping_quality < args.min_mapping_quality and after_read.is_supplementary) or \
            (after_read.has_tag('SA') and 'decoy' in after_read.get_tag('SA') or \
            (after_read.mapping_quality < args.min_mapping_quality and (not after_read.is_supplementary) and after_read.query_name in sup_reads)):
                reads_to_remove.add(after_read.query_name)


# add the reads without the ones we decided to remove
for after_read in after.fetch():
    if after_read.query_name not in reads_to_remove:
        output.write(after_read)

# add the reads from before
for before_read in before.fetch():
    if before_read.query_name in reads_to_remove:
        output.write(before_read)

# Close BAM files
before.close()
after.close()
output.close()