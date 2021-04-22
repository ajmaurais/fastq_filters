
import argparse
from multiprocessing import cpu_count

OPTIONS = ['query_offset', 'exclude_at_start', 'exclude_at_end']

FILTER_BASE_PARSER = argparse.ArgumentParser(add_help=False)

FILTER_BASE_PARSER.add_argument('-o', '--output_dir', default=None,
                                help='Directory to write output files to. Default is current workind directory.')
FILTER_BASE_PARSER.add_argument('--query_offset', default=0, type=int,
                                help='Offset at which expected_sequence should begin in query.')
FILTER_BASE_PARSER.add_argument('--exclude_at_start', default=0, type=int,
                                help='Length at beginning of alignment to ignore (length of diversity sequence)')
FILTER_BASE_PARSER.add_argument('--exclude_at_end', default=0, type=int,
                                help='Exclude any bases at the end? (e.g. because many truncated runs)')
FILTER_BASE_PARSER.add_argument('--nThread', default=int(cpu_count()), type=int,
                                help='Specify number of threads to use for parallel processing. '
                                     'Defaults to the number of logical cores in the system.')
FILTER_BASE_PARSER.add_argument('expected_sequence', help='Unmodified, WT sequence used to build library.')
FILTER_BASE_PARSER.add_argument('input_files', nargs='+', help='input file(s) to filter.')

