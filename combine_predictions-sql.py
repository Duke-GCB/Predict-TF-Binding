#!/usr/bin/env python
#
# combine_predictions.py
# Combines several prediction BED files into a single BED file
#
# When generating predictions from multiple cores, you'll have multiple prediction files that should be combined
# If the same region is found in multiple files, we only keep the highest score.

import argparse
import sqlite3

def read_bed_file(bed_file):
    for line in bed_file.readlines():
        yield line.strip().split()

def combine_predictions(bed_files):
    # Load all bed files into a database
    conn = sqlite3.connect(':memory:')
    conn.execute('''CREATE TABLE scores (chrom text, start int, end int, score real)''')
    for i, bed_file in enumerate(bed_files):
        for chrom, start, end, score in read_bed_file(bed_file):
            start, end, score = int(start), int(end), float(score)
            stmt = 'INSERT INTO scores VALUES (\'{}\',{},{},{})'.format(chrom, start, end, score)
            conn.execute(stmt)
        conn.commit()
    # Now use SQL to get the maximum score

    query = '''
    SELECT
        chrom,
        start,
        end,
        MAX(score) as score
    FROM
        scores
    GROUP BY
        chrom, start, end
    ORDER BY
        chrom ASC,
        start ASC
    '''
    for row in conn.execute(query):
        print ' '.join([str(x) for x in row])

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
