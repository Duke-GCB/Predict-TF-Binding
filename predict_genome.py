import itertools
import string
import subprocess
from svmutil import *

from Bio import SeqIO

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


def predict(matrix_file, model_file, output_file):
    """
    Lowest-level prediction function. Runs on files
    Runs svm-predict with the matrixfile and modelfile, writing results to resultsfile
    :param matrix_file:  Name of matrix file in libsvm format
    :param model_file:   Name of model file in libsvm format
    :param output_file: Name of file to store output of svm-predict
    :return: None
    """
    args = ["svm-predict", matrix_file, model_file, output_file]
    try:
        output, error = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        if len(error) > 0:
            raise Exception(error)
    except Exception as e:
        print "Error running svm-predict:", e
        raise e


def generate_matching_sequences(sequence, core, width):
    """
    Returns sub-sequences of width, that match the core in the middle.
    :param sequence: The the sequence to search, such as the whole sequence for a chromosome.
            Can be a string or a Bio.Seq
    :param core: The bases for which to search, in the center
    :param width: The desired sub-sequence width, e.g. 36
    :return: Generator, returning one sub-sequence per call
    """
    # Slide a window of width over the big sequence, and if the core is in the middle, return it
    # Subtract the window width to avoid indexing beyond the end of the sequence
    sequence_width = len(sequence)
    core_width = len(core)
    max_start = sequence_width - width
    # The core positions are calculated relative to the window (and not the overall sequence)
    # This works as long as both core and width are same parity
    core_start = (width - core_width) / 2
    for start in range(max_start):
        end = start + width
        window_sequence = sequence[start:end]
        window_core = window_sequence[core_start:core_start + core_width]
        if window_core == core:
            yield start, window_sequence

def read_genome_sequence(fasta_file, chrom):
    '''
    Reads a fasta file containing sequences for each chromosome in the genome, and returns the sequence for chrom
    :param fasta_file: File containing the genome sequences, such as hg19.fa
    :param chrom: chromosome name
    :return: The sequence, uppercased
    '''

    idx = SeqIO.index(fasta_file, 'fasta')
    record = idx[chrom]
    return record.seq.upper()


def load_model(model_file):
    return svm_load_model(model_file)


def predict_genome(genome_fasta_file, chrom, core, width, model_file, kmers):
    # 1. load the genome
    print "Loading {} from {}".format(chrom, genome_fasta_file)
    genome_sequence = read_genome_sequence(genome_fasta_file, chrom)
    # 2. Find the cores
    print "Generating matching sequences for core {}, width {}".format(core, width)
    matching_sequences = generate_matching_sequences(genome_sequence, core, width)
    # 3. Cores yield sequences
    # Currently evaluating libsvm bindings vs calling subprocess, so limiting to the first match
    for position, sequence in matching_sequences:
        print "Translating {} at position {} to SVR by kmers {}".format(str(sequence), position, kmers)
        # 4. Translate the sequences into SVR matrix by kmers
        features = svr_features_from_sequence(sequence, kmers)

        run_subprocess = True
        run_native = True

        if run_subprocess:
            matrixfile = 'matrix.txt'
            output_file = 'output.txt'
            print "writing features to matrixfile {}".format(matrixfile)
            write_svr_features([features], matrixfile)
            print "Predicting with results in {}".format(output_file)
            predict(matrixfile, model_file, output_file)

        if run_native:
            native_features = dict()
            # Build the dictionary that corresponds to the matrix file
            # Always has this 1:1 value at the beginning.
            native_features[1] = 1
            for i, feature in enumerate(features):
                native_features[i + 2] = feature['value']
            model = load_model(model_file) # Will move this out of the loop
            native_results = svm_predict([1], [native_features], model)
            print native_results

        break # Just do ONE.
        # TODO: Parse the prediction results, and join with the sequences / indices


if __name__ == '__main__':
    predict_genome('hg19.fa', 'chr16', 'GGAA', 36, 'ELK1_100nM_Bound_filtered_normalized_GGAA_1a2a3mer_format.model', [3])
