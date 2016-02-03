import itertools


def svr_features_from_sequence(seq, kmers):
    featureinfo = [['feature', 'position', 'featnum', 'featvalue'],
                   ['start', 'na', 1, 1]]  # header, and first feature which never changes
    featnum = 2  # the first feature is already defined as 1:1, so we start with 2
    for k in kmers:
        print "kmer: ", k
        bases = []  # creating empty lists needed later
        for n in itertools.product('ACGT', repeat=k): bases.append(''.join(n))  # all possible kmers with lenth k
        for n1 in range(len(seq) - (k - 1)):  # For each position in the sequence
            kmer = seq[n1:n1 + k]  # getting the actual kmer at this position
            for n2 in range(len(bases)):  # for every possible kmer of the size we want
                feature = bases[n2]  # getting the feature
                if feature == kmer:
                    featvalue = 1  # testing if the actual k-mer matches the feature, in which case the feature value is 1
                else:
                    featvalue = 0
                featureinfo.append([feature, n1, featnum, featvalue])  # adding the info about the feature to the list
                featnum += 1  # increasing the feature number by 1
    features = [0]  # starting a new list for building the matrix for SVR
    for x in range(1, len(featureinfo)):  # for every feature, in this list (skipping first item because header)
        features.append(
            str(featureinfo[x][2]) + ':' + str(featureinfo[x][3]))  # putting the feature values into the list
    return features


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
