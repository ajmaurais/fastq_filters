
import sys
import os
import argparse
from multiprocessing import cpu_count
import time

from .filter_generics import run_filter, FilterResult
from .base_parsers import FILTER_BASE_PARSER, MISMATCH_FILTER_BASE
from .base_parsers.mismatch_filter import OPTIONS as mismatch_filter_options
from .base_parsers.filter_base_parser import OPTIONS as filter_options
from .constants import QUERY_OFFSET, EXCLUDE_AT_END, EXCLUDE_AT_START, RANDOMIZED_LENGTH,\
        N_RANDOM_REGIONS, ADDITIONAL_MISMATCHES, QUERY_OFFSET


def mismatch_filter(query, expected_sequence,
                    query_offset=QUERY_OFFSET,
                    exclude_at_start=EXCLUDE_AT_START,
                    exclude_at_end=EXCLUDE_AT_END,
                    randomized_length=RANDOMIZED_LENGTH,
                    additional_mismatches=ADDITIONAL_MISMATCHES):
    '''
    Run base mismatch filter on query.

    Parameters
    ----------
    query: str
        Sequence to test.
    expected_sequence: str
        Non-ranzmized, base sequence.
    query_offset: int
        Offset at which expected_sequence should begin in query.
    exclude_at_start: int
        Length at beginning of alignment to ignore (length of diversity sequence)
    exclude_at_end: int
        Exclude any bases at the end? (e.g. because many truncated runs)
    randomized_length: int
        Length of ranzmized region.
    additional_mismatches: int
        Number of additional mismatches to allow.

    Returns
    -------
    results: FilterResult
    '''
    
    _query = query[query_offset:]
    begin_index = exclude_at_start
    end_index = len(expected_sequence) - exclude_at_end

    if len(_query) < end_index:
        return FilterResult(False, 'len', query)
    
    total_mismatches = 0
    consecutive_mismatches = 0
    longest_consecutive_mismtach = 0
    for i in range(begin_index, end_index):
        if _query[i] != expected_sequence[i]:
            total_mismatches += 1
            consecutive_mismatches += 1
            if consecutive_mismatches > longest_consecutive_mismtach:
                longest_consecutive_mismtach = consecutive_mismatches
        else:
            consecutive_mismatches = 0

    total_mismatches -= min((randomized_length, longest_consecutive_mismtach))
    if total_mismatches > additional_mismatches:
        return FilterResult(False, 'mismatch', query)    
    return FilterResult(True, None, query)


def sequence_generator(fname):
    with open(fname, 'r') as inF:
        for line in inF:
            yield line.strip()


def run_mismatch_filter(input_file, expected_sequence,
                        output_dir_path=None,
                        nThread=1, chunk_size=int(1e6), **kwargs):

    # initialize output file names
    fastq_base = os.path.splitext(os.path.basename(input_file))[0]
    output_dir = os.path.abspath(os.path.dirname(input_file)) if output_dir_path is None else output_dir_path
    output_file_name = '{}/{}.mismatchFiltered.txt'.format(output_dir, fastq_base)
    failed_file_name = '{}/{}.mismatchFailed.txt'.format(output_dir, fastq_base)

    # start timer
    start_time = time.time()
    sys.stdout.write('Working on {}...\n'.format(input_file))
    sys.stdout.write('\tStarted filtering reads at {}\n'.format(time.strftime('%H:%M:%S',
                                                                time.localtime())))
    # Get entry generator
    seq_gen = sequence_generator(input_file)

    filter_results = run_filter(output_file_name, failed_file_name,
                                seq_gen, mismatch_filter, expected_sequence,
                                verbose=True, use_pool=False, **kwargs)
    # Print summary statistics to stdout
    out = sys.stdout
    total = filter_results.total()

    out.write('\n\tQ score filtering results for {}\n'.format(fastq_base))
    out.write('\t{:,} reads processed\n'.format(total))
    out.write('\t{:,} were the incorrect length ({:.2%})\n'.format(filter_results.results_count['len'],
                                                                  filter_results.results_count['len'] / total))
    out.write('\t{:,} had too many base mismatches ({:.2%})\n'.format(filter_results.results_count['mismatch'],
                                                                      filter_results.results_count['mismatch'] / total))
    out.write('\tFiltered results written to {}\n'.format(output_file_name))

    # print elapsed time
    end_time = time.time()
    out.write('Total run time = {:.0f} seconds ({:.1f} minutes)\n'.format(end_time - start_time,
                                                                          (end_time - start_time) / 60))


def main():
    parser = argparse.ArgumentParser(description='Apply mismatch filter to txt files.',
                                     parents=[FILTER_BASE_PARSER, MISMATCH_FILTER_BASE])
    args = parser.parse_args()

    for fname in args.input_files:
        run_mismatch_filter(fname, args.expected_sequence,
                            nThread=args.nThread,
                            chunk_size=1e6,
                            update_period=int(1e6),
                            **{x: getattr(args, x) for x in mismatch_filter_options + filter_options})


if __name__ == '__main__':
    main()


