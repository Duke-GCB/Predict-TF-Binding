import itertools

NUCLEOTIDES='ACGT'

def svr_features_from_sequence(seq, kmers):
    """
    Transforms a sequence and list of kmer values into a list of dictionaries
    that be converted easily into an SVR matrix.

    kmers is a list of lengths. For each length, the function will enumerate all possible
    nucleotide combinations at that length ('AA','AC','AG',...'TT')

    The function returns a list of all possible positions in the sequence x all possible features
    and indicates 1 if the feature matches that position in the sequence, or 0 if it does not.

    For example, for the input sequence 'ACAGTC' and a kmer value of [2,3], the function
     produces the following

    [{'position': 0, 'value': 0, 'feature': 'AA'}, # 'AA' does not match at position 0
     {'position': 0, 'value': 1, 'feature': 'AC'}, # 'AC' matches at position 0 in 'ACAGTC'
     {'position': 0, 'value': 0, 'feature': 'AG'},
     ...
     {'position': 1, 'value': 1, 'feature': 'CA'}, # 'CA' matches at position 1 in 'ACAGTC'
     ...
     {'position': 0, 'value': 0, 'feature': 'AAT'},
     {'position': 0, 'value': 1, 'feature': 'ACA'},
     {'position': 0, 'value': 0, 'feature': 'ACC'},

    :param seq: A sequence of nucleotides to expand
    :param kmers: list of integers (e.g. [1,2,3])
    :return: a list of dictionaries, containing position, featvalue, and feature
    """

    svr_features = []
    for k in kmers:
        # Generate all possible combinations of length k (e.g ['AAA', 'AAC', ...  'TTG', 'TTT']
        features = [''.join(x) for x in itertools.product(NUCLEOTIDES, repeat=k)]
        # Check each position in the sequence for a match
        n_sub_seqs = len(seq) - (k - 1) # If seq length is 36 and k is 3, there are 34 positions
        for position in range(n_sub_seqs):
            sub_seq = seq[position:position + k] # the sub-sequence with length k
            for feature in features:
                if feature == sub_seq:
                    value = 1
                else:
                    value = 0
                info = {'feature': feature, 'position': position, 'value' : value}
                svr_features.append(info)
    return svr_features


def svr_matrix_from_seqlist(seqlist, kmers):
    svrmatrix = []
    for seq in seqlist:
        svrmatrix.append(svr_features_from_sequence(seq, kmers))
    return svrmatrix


def get_core(seq):
    return seq[16:20]


kmers = [3]
seqlist = ['GCCCGCAGAGCGGAAGGCGGGATGGCTGGGGGCGGG',
           'GGGCTCAGCGCCGACTGCGCGCCTCTGCCCGCGAAA',
           'CGCCATAGCGACGGCGCCGCAATTTAGGAGCGTGCT',
           'TCCCGCCTTCCGCTCTGCGGGCGGCAGCCGGGCTGG']

matrix = svr_matrix_from_seqlist(seqlist, [3])
for seq in matrix:
    for x in seq:
        print x
