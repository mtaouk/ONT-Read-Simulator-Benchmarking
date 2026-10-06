#!/usr/bin/env python3

import argparse
import csv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('tool', help='Tool name in the TSV, e.g. pbsim3_hp_5')
    parser.add_argument('tsv', help='Homopolymer results TSV containing real and tool results')
    parser.add_argument('counts', help='Comma-separated reference counts for lengths 3, 4, 5, ...')
    args = parser.parse_args()

    try:
        counts = [int(value) for value in args.counts.split(',')]
    except ValueError:
        parser.error('Reference counts must be comma-separated integers')
    if any(count < 0 for count in counts) or sum(counts) == 0:
        parser.error('Reference counts must be nonnegative with a positive total')

    with open(args.tsv, newline='') as tsv:
        error_rates = {
            (row['tool'], int(row['rlen'])): 1 - float(row['ATGC'])
            for row in csv.DictReader(tsv, delimiter='\t')
        }

    weighted_difference = 0.0
    for hp_len, count in enumerate(counts, start=3):
        for tool in ('real', args.tool):
            if (tool, hp_len) not in error_rates:
                parser.error(f'Missing result for {tool}, homopolymer length {hp_len}')
        weighted_difference += count * abs(
            error_rates[args.tool, hp_len] - error_rates['real', hp_len]
        )

    print(weighted_difference / sum(counts))


if __name__ == '__main__':
    main()
