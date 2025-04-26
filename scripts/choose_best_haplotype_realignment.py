#!/usr/bin/env python3
import argparse
import itertools
import pysam
import numpy as np
import tqdm.auto as tqdm
import logging
import pandas as pd

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

    
def _poor_mq_has_sa(mappings: pd.DataFrame, mq_threshold: int) -> bool:
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
    if len(mappings) > 1:
        has_sa = True
    else:
        has_sa = False
    if not has_sa:
        return False
    
    return np.any(~mappings['is_supplementary'] & mappings['mapq'] < mq_threshold)
    
def _poor_mq(mappings: list) -> bool:
    return np.any(mappings['poor_sa'])


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
    return max([x.qend for x in mappings]) - min([x.qstart for x in mappings])

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
    first_poor_mq =  _poor_mq_has_sa(first, mq_threshold) or _poor_mq(first)
    second_poor_mq =  _poor_mq_has_sa(second, mq_threshold) or _poor_mq(second)
    
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
    
def find_best_mapping_index(mappings: list, tie_breaker_idx: int = None) -> int:
    '''Find the best mapping index

    Parameters
    ----------
    mappings : list
        List of dataframes representing alignments of the same read
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
        comp = compare_read_mappings(mappings[best_index], mappings[i], best_index, i, tie_breaker_idx=tie_breaker_idx)
        if comp == -1:
            best_index = i
    return best_index

def _chunk_to_df(chunk: iter, mq_threshold: int):
    '''Convert a chunk of mappings to a dataframe

    Parameters
    ----------
    chunk : iter
        Chunk of mappings

    Returns
    -------
    pd.DataFrame
        Dataframe of mappings
    '''
    qvals = [ (x.query_name, x.reference_id, x.query_alignment_start, 
               x.query_alignment_end, x.is_supplementary, x.mapping_quality, 
               x.get_tag("SA") if x.has_tag("SA") else "" ) for x in chunk ]
 
    df = pd.DataFrame(qvals)
 
    df[1] = df[1].astype('int16')
    df[5] = df[5].astype('int8')
    df[2] = df[2].astype('int16')
    df[3] = df[3].astype('int16')
    df[6] = ((df[4] & df[5] < mq_threshold) | (df[6].str.contains('decoy'))).astype(bool)
    df[0] = df[0].astype('string[pyarrow]')
    df.columns = ['read_name', 'chrom','qstart', 'qend','is_supplementary','mapq','poor_sa']
    return df

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return itertools.zip_longest(fillvalue=fillvalue, *args)

def get_alignment_df( input_file: str, mq_threshold: int):
    '''Get the alignment dataframe from the input file

    Parameters
    ----------
    input_file : str
        Input file name

    Returns
    -------
    pd.DataFrame
        Dataframe of mappings
    '''
    dfs = []
    with pysam.AlignmentFile(input_file, "rb") as input:
        for chunk in tqdm.tqdm(grouper(input, 500000)):
            df = _chunk_to_df(chunk, mq_threshold)
            dfs.append(df)
    return pd.concat(dfs, ignore_index=True).sort_index()

def run():
    args = parse_args()
    logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
    logger = logging.getLogger(__name__)    
    logger.info("Merging alignments: start")
    logger.info(f"Alignment sources: {args.alignment_sources}")
    logger.info(f"Output: {args.output}")
    logger.info(f"Min mapping quality: {args.min_mapping_quality}")
    logger.info(f"Tie breaker index: {args.tie_breaker_idx}")
    logger.info("Reading alignment files and creating mapping keys:start")

    mappings = [ get_alignment_df(x).groupby("read_name") for x in args.alignment_sources ]
    logger.info("Reading alignment files and creating mapping keys:done")
    best_mappings = {}
    logger.info("Choosing best alignment per read: start")    
    for k in mappings[0].groups.keys():
        best_mappings[k] = find_best_mapping_index([x.get_group(k) for x in mappings], args.tie_breaker_idx)
    logger.info("Choosing best alignment per read: end")    


    # Create a new BAM file for output
    counters = np.zeros(len(args.alignment_sources))
    logger.info("Reading from the sources and writing into the output:start")
    with pysam.AlignmentFile(args.output, "wb", template=pysam.AlignmentFile(args.alignment_sources[0])) as output:
        for i, r in enumerate(args.alignment_sources):
            logger.info(f"Reading alignments from {r}")
            input = pysam.AlignmentFile(r, "rb")
            for rec in tqdm.tqdm(input):
                if best_mappings[rec.query_name] == i:
                    if args.overwrite_mapq is not None and args.overwrite_mapq[0] == i:
                        if rec.mapping_quality >= args.overwrite_mapq[1] and rec.mapping_quality < args.overwrite_mapq[2]:
                            rec.mapping_quality = args.overwrite_mapq[2]
                    output.write(rec)
                    counters[i] += 1
            logger.info(f"Wrote {counters[i]} alignments from {r}")
    output.close()

if __name__ == "__main__":
    run()                       