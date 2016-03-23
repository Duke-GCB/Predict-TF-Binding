#!/usr/bin/env python
#
# resize_ranges
# Resizes ranges in a  BED file around a center point
#
# When generating predictions from multiple cores, you'll have multiple prediction files that should be combined
# If the same region is found in multiple files, we only keep the highest score.

import argparse
import sqlite3

def read_bed_file(bed_file):
    for line in bed_file.readlines():
        yield line.strip().split()

def resize_range(start, end, new_width):
    start, end = int(start), int(end)
    original_range = end-start
    margin = (original_range - new_width) / 2
    start = start + margin
    end = end - margin
    return start, end

def resize_ranges(input_file, width, output_file_name):
    with open(output_file_name, 'w') as f:
        for chr, start, end, score in read_bed_file(input_file):
            start, end = resize_range(start, end, width)
            print >> f, chr, start, end, score

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Resize ranges in a BED file')
    parser.add_argument('input_file',
                        metavar='input_file',
                        type=file,
                        help='Bed file to consider',)
    parser.add_argument('output_file',
                        metavar='output_file',
                        help='Output bed file',)
    parser.add_argument('--width',
                        metavar='width',
                        type=int,
                        default=20,
                        help='Desired width of the range',)

    args = parser.parse_args()
    resize_ranges(args.input_file, args.width, args.output_file)
