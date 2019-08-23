#!/usr/bin/env python

import sys
import argparse
import csv


def change_precision(input, output, precision, delimiter, source_index=3):
    """
    Changes the precision of the value at source_index to precision
    value of one column by a factor
    :param input: An input stream or open file
    :param output: An output stream
    :param precision: Number of decimal places to use
    :param delimiter: Separator for the file, e.g. tab, space, or comma.
    :param source_index: Column index containing the source value
    :return: None
    """
    reader = csv.reader(input, delimiter=delimiter)
    writer = csv.writer(output, delimiter=delimiter)
    for row in reader:
        orig_prediction = float(row[source_index])
        row[source_index] = '{:.{prec}f}'.format(orig_prediction, prec=precision)
        writer.writerow(row)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Change precision of a float column in a BED file')
    parser.add_argument('inputfile', type=argparse.FileType('r'))
    parser.add_argument('precision', type=int)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--spaces', action='store_const', const=' ', dest='delimiter')
    group.add_argument('--tabs', action='store_const', const='\t', dest='delimiter')
    group.add_argument('--commas', action='store_const', const=',', dest='delimiter')
    parser.add_argument_group()
    args = parser.parse_args()
    change_precision(args.inputfile, sys.stdout, args.precision, args.delimiter)
