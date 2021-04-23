
import os
import re

from .filter_generics import print_alignment
from .base_parsers.filter_base_parser import OPTIONS as filter_options

def write_allignment(fname, expected_sequence, query_offset=0):
    ofname = '{}.alignment.txt'.format(os.path.splitext(fname)[0])

    with open(fname, 'r') as inF:
            lines = [x.strip() for x in inF.readlines()]

    with open(ofname, 'w') as outF:
        for line in lines:
                _query = line[query_offset:]
            begin_index = exclude_at_start
            end_index = len(expected_sequence) - exclude_at_end
            print_alignment(expected_sequence, _query, out=outF)


def main():
    parser = argparse.ArgumentParser(description='Apply Q filter to fastq file(s).',
                                 parents=[FILTER_BASE_PARSER])

    parser.add_argument('--query_offset', default=0, type=int,
                    help='Offset at which expected_sequence should begin in query.')

    parser.add_argument('expected_sequence', help='Expected sequence. Must match the regex "{}"'.format(SEQUENCE_RE))

    FILTER_BASE_PARSER.add_argument('input_files', nargs='+', help='input file(s) to filter.')

    args = parser.parse_args()

    if re.match(args.expected_sequence):
        _expected_sequence = args.expected_sequence
    else:
        raise RuntimeError('{} is an invalid sequence!'.format(SEQUENCE_RE))

    for fname in args.input_files:
        write_allignment(fname, _expected_sequence, args.query_offset)
