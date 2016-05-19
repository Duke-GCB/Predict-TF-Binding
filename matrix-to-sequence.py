#!/usr/bin/env python

# Script to convert a matrix back to a sequence for investigation
# Currently hard-coded to assume 20 1-mers at the beginning of the matrix

import fileinput

MAP= {
  '1000': 'A',
  '0100': 'C',
  '0010': 'G',
  '0001': 'T',
}

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

for line in fileinput.input():
  sequence = []
  line = line.strip()
  score = float(line.split()[0])
  terms = line.split()[1:81]
  for c in chunks(terms, 4):
    code = ''.join([x.split(':')[1] for x in c])
    sequence.append(MAP[code])
  print ''.join(sequence), score

