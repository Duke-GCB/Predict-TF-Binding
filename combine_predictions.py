#!/usr/bin/env python
#
# combine_predictions.py
# Combines several prediction BED files into a single BED file
#
# When generating predictions from multiple cores, you'll have multiple prediction files that should be combined
# If the same region is found in multiple files, we only keep the highest score.

import argparse

def read_bed_file(bed_file):
    for line in bed_file.readlines():
        yield line.strip().split()

def combine_predictions(bed_files):
    # (chr1, 1, 2)
    max_values = {}
    for bed_file in bed_files:
        for chrom, start, end, score in read_bed_file(bed_file):
            start, end, score = int(start), int(end), float(score)
            k = (chrom, start, end)
            if k in max_values.keys():
                previous_score = max_values[k]
                max_values[k] = max(score, previous_score)
            else:
                max_values[k] = score
    for k, v in sorted(max_values.items()):
        print k[0], k[1], k[2], v

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Combines several prediction BED files into a single BED file, using'
                                                   'the highest score in the case of duplicate regions')
    parser.add_argument('bedfiles',
                        metavar='bed_file',
                        type=file,
                        help='Bed files to consider',
                        nargs='+',)
    args = parser.parse_args()
    combine_predictions(args.bedfiles)
