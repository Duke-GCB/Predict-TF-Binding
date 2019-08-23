#!/usr/bin/env python

import argparse
import csv
import sys


def filter_scores(input, output, delimiter, threshhold=0.0, source_index=3):
    """
    Filters a predictions bed file by returning only rows where the score is
    above the threshold
    :param input: An input stream or open file
    :param output: An output stream
    :param delimiter: input and output file delimiter
    :param threshhold: minimum value for inclusion
    :param source_index: Column index containing the source value
    :return:
    """
    reader = csv.reader(input, delimiter=delimiter)
    writer = csv.writer(output, delimiter=delimiter)
    for row in reader:
        # Adds a score column by multiplying the value of an existing column by a factor
        # http://genome.ucsc.edu/FAQ/FAQformat.html#format1
        if float(row[source_index]) > threshhold:
            writer.writerow(row)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfile', type=argparse.FileType('r'))
    parser.add_argument('threshhold', type=float, default=0.0)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--spaces', action='store_const', const=' ', dest='delimiter')
    group.add_argument('--tabs', action='store_const', const='\t', dest='delimiter')
    group.add_argument('--commas', action='store_const', const=',', dest='delimiter')
    parser.add_argument_group()
    args = parser.parse_args()
    filter_scores(args.inputfile, sys.stdout, args.delimiter, args.threshhold)
