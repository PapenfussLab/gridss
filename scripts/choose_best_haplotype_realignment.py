#!/usr/bin/env python3
import argparse
import itertools
import pysam
import numpy as np
import tqdm.auto as tqdm
import logging

def comma_sep(x):
    xsp = x.split(',')
    return (xsp[0],xsp[1],xsp[2])
def parse_args():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Merge haplotype alignments from multiple sources')
    parser.add_argument('--alignment_source', required=True, type=str, action='append', help='Alignment sources (specify multiple times)')
    parser.add_argument('--output', required=True, type=str,help='The output realigned file')
    parser.add_argument('--min_mapping_quality', type=int, default=60, help='Supplementary alignments with mapping quality below this threshold will be considered poor')
    parser.add_argument('--tie_breaker_idx', type=int, default=None, help='Index of the alignment source to prioritize in case of a tie (e.g. giraffe)')
    parser.add_argument('--overwrite_mapq', type=comma_sep, help="Comma-separated tuple of three values: index of source to overwrite, min_mapq to change, value_to_assign. Useful for overwriting low-scale mapping qualities of giraffe")
    args = parser.parse_args()
    args.alignment_sources = args.alignment_source
    return args

class MappingKey:
     # keep read name, chromosome, read start, read end, mapping quality and if the alignment is supplementary
    def __init__(self, rec:pysam.AlignedSegment):
        self.name = rec.query_name
        self.chrom = rec.reference_name
        self.start = rec.query_alignment_start
        self.end = rec.query_alignment_end
        self.mapq = rec.mapping_quality
        self.is_supplementary = rec.is_supplementary
        self.sa_tag = None
        if rec.has_tag('SA'):
            self.sa_tag = rec.get_tag('SA')

    def poor_sa(self, mq_threshold) -> bool:
        if self.sa_tag is not None and 'decoy' in self.sa_tag:
            return True
        if self.is_supplementary and self.mapq < mq_threshold:
            return True
        return False 
    
def _poor_mq_has_sa(mappings: list, mq_threshold: int) -> bool:
    '''Check if there is a poor mapping quality read with supplementary alignment

    Parameters
    ----------
    mappings : list
        List of MappingKey objects

    Returns
    -------
    bool
        True if there is a poor mapping quality read with supplementary alignment, False otherwise
    '''
    has_sa = False
    for mapping in mappings:
        if mapping.is_supplementary:
            has_sa = True
            break
    poor_mq = False
    if has_sa:
        for mapping in mappings:
            if not mapping.is_supplementary and mapping.mapq < mq_threshold:
                poor_mq = True
                break
            
    return poor_mq

def _q_span(mappings: list) -> int:
    '''Calculate the query span of the mappings

    Parameters
    ----------
    mappings : list
        List of MappingKey objects

    Returns
    -------
    int
        Query span
    '''
    return max([x.end for x in mappings]) - min([x.start for x in mappings])

def compare_read_mappings(first: list, second: list, mq_threshold: int, idx1:int, idx2: int, tie_breaker_idx: int) -> int:
    '''Compare two lists of MappingKey objects. The comparison is based on the following criteria:
    1. If one of the alignments is from a source with poor mapping quality, we choose the one from the source with better quality
    2. If both alignments are from sources with poor mapping quality, we consider them equal
    3. If one alignment is not split into two parts and the other is, we choose the one that is not split
    4. If both alignments are split into two parts, we choose the one with the larger query span
    5. If the query span is equal, we choose the one from the source specified as a tie breaker

    Parameters
    ----------
    first : list
        List of MappingKey objects
    second : list
        List of MappingKey objects
    mq_threshold : int
        Mapping quality threshold
    idx1 : int
        Sourse of the first alignment
    idx2 : int
        Sourse of the second alignment
    tie_breaker_idx : int
        If one of the alignments is from `tie_breakder_idx` source (e.g. giraffe) and they are of equal quality - we choose the one from `tie_breaker_idx`

    Returns
    -------
    int
        1 if first > second, -1 if second > first, 0 if equal
    '''
    first_poor_mq = _poor_mq_has_sa(first, mq_threshold) or np.any([x.poor_sa(mq_threshold) for x in first])
    second_poor_mq = _poor_mq_has_sa(second, mq_threshold) or np.any([x.poor_sa(mq_threshold) for x in second])
    if not first_poor_mq and second_poor_mq:
        return 1
    if first_poor_mq and not second_poor_mq:
        return -1   
    if first_poor_mq and second_poor_mq:
        return 0
    
    first_q_span = _q_span(first)
    second_q_span = _q_span(second)

    # prioritize cases when the mapping is not split into two parts
    if len(first) > 1 and len(second) == 1:
        if first_q_span - 10 <= second_q_span:
            return -1
    if len(first) == 1 and len(second) > 1:
        if second_q_span - 10 <= first_q_span:
            return 1
    
    if first_q_span > second_q_span:
        return 1
    if first_q_span < second_q_span:
        return -1
    if idx1 == tie_breaker_idx:
        return 1
    if idx2 == tie_breaker_idx:
        return -1
    return 0
    
def find_best_mapping_index(mappings: list, mq_threshold: int, tie_breaker_idx: int = None) -> int:
    '''Find the best mapping index

    Parameters
    ----------
    mappings : list
        List of MappingKey objects
    mq_threshold : int
        Mapping quality threshold
    tie_breaker_idx : int
        Provided if one of the alignments is from `tie_breakder_idx` source (e.g. giraffe)
    Returns
    -------
    int
        Best mapping index
    
    See also
    --------
    compare_read_mappings
    '''
    best_index = 0
    for i in range(1, len(mappings)):
        comp = compare_read_mappings(mappings[best_index], mappings[i], mq_threshold, best_index, i, tie_breaker_idx=tie_breaker_idx)
        if comp == -1:
            best_index = i
    return best_index

def validate_sort_order(alignment_sources: list):
    '''Validate that the sort order of the alignment sources is by queryname

    Parameters
    ----------
    alignment_sources : list
        List of alignment sources

    Raises
    ------
    ValueError
        If the sort order of the alignment sources is not the same
    '''
    sort_orders = [pysam.AlignmentFile(x, "rb").header.get('HD', {}).get('SO') for x in alignment_sources]
    for i,so in enumerate(sort_orders):
        if so != 'queryname':
            raise ValueError(f"Alignment source {alignment_sources[i]} is not sorted by queryname")

def run():
    args = parse_args()
    logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
    logger = logging.getLogger(__name__)    
    logger.info("Merging alignments: start")
    logger.info(f"Alignment sources: {args.alignment_sources}")
    logger.info(f"Output: {args.output}")
    logger.info(f"Min mapping quality: {args.min_mapping_quality}")
    logger.info(f"Tie breaker index: {args.tie_breaker_idx}")
    validate_sort_order(args.alignment_sources)
    mappings = zip(*(itertools.groupby(map(MappingKey, pysam.AlignmentFile(x, "rb")), lambda x: x.name) for x in args.alignment_sources))
    logger.info("Choosing best alignment per read: start")    
    counters = np.zeros(len(args.alignment_sources))
    with pysam.AlignmentFile(args.output, "wb", template=pysam.AlignmentFile(args.alignment_sources[0])) as output:
        for alns in tqdm.tqdm(mappings):
            assert [x[0] for x in alns] == [alns[0][0]] * 2, "Read names are not equal"
            alignments = [list(x[1]) for x in alns]
            best_mapping_idx = find_best_mapping_index(alignments, args.min_mapping_quality, args.tie_breaker_idx)
            best_alignments = alignments[best_mapping_idx]
            for rec in best_alignments:
                if args.overwrite_mapq is not None and args.overwrite_mapq[0] == best_mapping_idx:
                    if rec.mapping_quality >= args.overwrite_mapq[1] and rec.mapping_quality < args.overwrite_mapq[2]:
                        rec.mapping_quality = args.overwrite_mapq[2]
                    output.write(rec)
            counters[best_mapping_idx] += 1
    for i in range(len(args.alignment_sources)):
        logger.info(f"Wrote {counters[i]} alignments from {args.alignment_sources[i]} to {args.output}")
    output.close()
    logger.info("Choosing best alignment per read: start")    

if __name__ == "__main__":
    run()                       