import itertools

NUCLEOTIDES='ACGT'

def svr_features_from_sequence(seq, kmers):
    featureinfo = [['feature', 'position', 'featnum', 'featvalue'],
                   ['start', 'na', 1, 1]]  # header, and first feature which never changes
    for k in kmers:
        print "kmer: ", k
        # Generate all possible combinations of length k (e.g ['AAA', 'AAC', ...  'TTG', 'TTT']
        features = [''.join(x) for x in itertools.product(NUCLEOTIDES, repeat=k)]
        # Check each position in the sequence for a match
        n_sub_seqs = len(seq) - (k - 1) # If seq length is 36 and k is 3, there are 34 positions
        for position in range(n_sub_seqs):
            sub_seq = seq[position:position + k] # the sub-sequence with length k
            for feature_index, feature in enumerate(features):
                if feature == sub_seq:
                    feature_value = 1
                else:
                    feature_value = 0
                info = {'feature': feature, 'position': position, 'featnum': feature_index + 2, 'featvalue' : feature_value}
                featureinfo.append(info)
    return featureinfo


def svr_matrix_from_seqlist(seqlist, kmers):
    ### Generating matrix file for each sequence in this set
    svrmatrix = []
    for seq in seqlist:
        svrmatrix.append(svr_features_from_sequence(seq, kmers))
    return svrmatrix


def get_core(seq):
    return seq[16:20]


# questions: how does the kmers list get interpreted

kmers = [3]
seqlist = ['GCCCGCAGAGCGGAAGGCGGGATGGCTGGGGGCGGG',
           'GGGCTCAGCGCCGACTGCGCGCCTCTGCCCGCGAAA',
           'CGCCATAGCGACGGCGCCGCAATTTAGGAGCGTGCT',
           'TCCCGCCTTCCGCTCTGCGGGCGGCAGCCGGGCTGG']

print "Cores:", ','.join([get_core(seq) for seq in seqlist])

modelfile = 'ELK1_100nM_Bound_filtered_normalized_GGAA_1a2a3mer_format.model'
matrix = svr_matrix_from_seqlist(seqlist, [3])
print matrix
