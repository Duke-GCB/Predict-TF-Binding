import itertools
import string
import subprocess

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

def write_svr_features(all_svr_features, matrixfile):
    """
    Writes a list of lists of feature dictionaries to a file
    :param all_svr_features: A list of lists of svr_feature dictionaries. One entry in the master list per sequence
    :param matrixfile: file name to write to
    :return: None
    """
    with open(matrixfile, 'w') as f:
        for svr_features in all_svr_features:
            print >> f, "0\t1:1\t", # Required header for matrix file
            # enumerate yields tuples with index (0, int) and the item (1, dict)
            colon_separated = map(lambda x:  '{}:{}'.format(x[0] + 2,x[1]['value']), enumerate(svr_features))
            print >> f, '\t'.join(colon_separated)


def load_sequences(sequence_file):
    """
    Loads sequences from a text file, one per line
    :param sequence_file: a text file, with one sequence per line
    :return: list of sequence strings
    """
    with open(sequence_file,'r') as f:
        sequences = map(string.strip, f)
    return sequences


def predict(matrixfile, modelfile, outputfile):
    """
    Lowest-level prediction function. Runs on files
    Runs svm-predict with the matrixfile and modelfile, writing results to resultsfile
    :param matrixfile:  Name of matrix file in libsvm format
    :param modelfile:   Name of model file in libsvm format
    :param outputfile: Name of file to store output of svm-predict
    :return: None
    """
    args = ["svm-predict", matrixfile, modelfile, outputfile]
    try:
        output, error = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        if len(error) > 0:
            raise Exception(error)
    except Exception as e:
        print "Error running svm-predict:", e
        raise e


if __name__ == '__main__':
    sequences = ['GCCCGCAGAGCGGAAGGCGGGATGGCTGGGGGCGGG',
                 'GGGCTCAGCGCCGACTGCGCGCCTCTGCCCGCGAAA',
                 'CGCCATAGCGACGGCGCCGCAATTTAGGAGCGTGCT',
                 'CAGGCTTTGGGAGCCAGCGGGGCGGGAGCGGCGAAG']
    all_features = map(lambda seq: svr_features_from_sequence(seq, [3]), sequences)
    write_svr_features(all_features, 'svr_matrix.txt')
