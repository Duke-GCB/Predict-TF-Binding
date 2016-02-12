import itertools
import string
import os
import sys
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
    str_seq = str(seq) # If seq is a Bio.Seq, it's faster to check it as a string
    svr_features = []
    for k in kmers:
        # Generate all possible combinations of length k (e.g ['AAA', 'AAC', ...  'TTG', 'TTT']
        features = [''.join(x) for x in itertools.product(NUCLEOTIDES, repeat=k)]
        # Check each position in the sequence for a match
        n_sub_seqs = len(str_seq) - (k - 1) # If seq length is 36 and k is 3, there are 34 positions
        for position in range(n_sub_seqs):
            sub_seq = str_seq[position:position + k] # the sub-sequence with length k
            # start with a template list. All zero values at the current position
            exploded = [{'feature': feature, 'position': position, 'value': 0} for feature in features]
            # Determine the index of the current sub seq
            try:
                feature_index = features.index(sub_seq)
                exploded[feature_index]['value'] = 1
            except ValueError:
                print "Warning: sub-sequence '{}' not found in features".format(sub_seq)
            svr_features.extend(exploded)
    return svr_features


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

def read_genome_idx(fasta_file):
    """
    Reads a fasta-formatted genome and returns a Bio.SeqIO.index
    :param fasta_file: File containing the genome sequences, such as hg19.fa
    :return: A Bio.SeqIO index
    """
    return SeqIO.index(fasta_file, 'fasta')


def get_chrom_sequence(idx, chrom):
    """
    Gets the sequence for a chromosome
    :param idx: The BioPython sequence, indexed by chrom
    :param chrom: chromosome name
    :return: The sequence, uppercased
    """
    record = idx[chrom]
    return record.seq.upper()


def load_model(model_file):
    """
    Loads a svm model from a file and computes its size
    :param model_file: The file name of the model to load
    :return: A dictionary with keys model, file, and size
    """
    model = svm_load_model(model_file)
    size = len(model.get_SV()[0])
    model_dict = {'model': model, 'size': size}
    return model_dict


def predict(features, model, const_intercept=False):
    """
    Run prediction using svm_predict.
    :param features: List of features, produced by svr_features_from_sequence
    :param model: A loaded svm model (from load_model)
    :param const_intercept: if true, add a 1:1 term at the beginning of the matrix. Must match model's term
    :return: triple of predictions, accuracy, and values from svm_predict.
    """
    svm_matrix = dict()
    # Build the dictionary that corresponds to the matrix file
    offset = 1 # svm_matrix is a dictionary of index to value, starting at 1
    if const_intercept:
        svm_matrix[offset] = 1
        offset += 1
    for i, feature in enumerate(features):
        svm_matrix[i + offset] = feature['value']
    # svm_predict defines an info function that prints results to STDOUT
    # So we suppress this by temporarily assigning sys.stdout to os.devnull
    old_stdout = sys.stdout
    devnull = open(os.devnull, 'w')
    sys.stdout = devnull
    predictions = svm_predict([1], [svm_matrix], model)
    # Restore sys.stdout
    sys.stdout = old_stdout
    return predictions


def print_bed(file_handle, chrom, position, width, score):
    """
    Prints an annotation in BED format (space-separated) to the provided file handle, could be sys.stdout or an open file
    :param file_handle: destination file handle for line printing
    :param chrom: Chromosome, e.g. 'chr16'
    :param position: Integer position of the annotation
    :param width: Integer width of the annotation
    :param score: floating-point value of the annotation
    :return: None
    """
    print >> file_handle, chrom, position, position + width, score

def predict_genome(genome_fasta_file, core, width, model_file, kmers, const_intercept, output_file):
    """
    Generate predictions on the provided genome fasta file.
    Predictions will only be generated on sequences of width 'width', matching the nucleotides of 'core' in the center

    :param genome_fasta_file: File name of fasta-formatted genome, e.g. hg19.fa or hg38.fa
    :param core: sequence of nucleotides to find when generating predictions
    :param width: width, in bases, of the window on which to generate predictions
    :param model_file: Name of the svm model file to load
    :param kmers: List of integers (e.g. [1,2,3]) for base combination in prediction generation. Must match model generation parameters
    :param const_intercept: true or false - whether or not to add a constant term to the matrix generation. Must match model generation
    :param output_file: Output file to write, in bed format
    :return: None
    """
    # 1. load the genome
    print 'Loading genome', genome_fasta_file
    idx = read_genome_idx(genome_fasta_file)

    # 2. Load model
    print 'Loading model', model_file
    model_dict = load_model(model_file)

    # 3. Iterate over all chromosomes in genome
    with open(output_file, 'w') as output:
        for chrom in idx:
            print 'Predicting on', chrom
            # Run prediction for the chrom
            for position, sequence, score in predict_chrom(idx, chrom, core, width, model_dict, kmers, const_intercept):
                print_bed(sys.stdout, chrom, position, width, score)
                print_bed(output, chrom, position, width, score)


def predict_chrom(sequence_idx, chrom, core, width, model_dict, kmers, const_intercept):
    """
    Generates predictions for a single chromosome
    :param sequence_idx: The indexed SeqIO object
    :param chrom: The name of a chromosome (must be the id in the index)
    :param core: sequence of nucleotides to find when generating predictions
    :param width: width, in bases, of the window on which to generate predictions
    :param model_file: Name of the svm model file to load
    :param kmers: List of integers (e.g. [1,2,3]) for base combination in prediction generation. Must match model generation parameters
    :param const_intercept: true or false - whether or not to add a constant term to the matrix generation. Must match model generation
    :return: A generator, yielding the start position, sequence, and prediction score
    """

    chrom_sequence = get_chrom_sequence(sequence_idx, chrom)
    print "Generating matching sequences for core {}, width {}".format(core, width)
    matching_sequences = generate_matching_sequences(chrom_sequence, core, width)
    # 3. Cores yield sequences

    for position, sequence in matching_sequences:
        # 4. Translate the sequences into SVR matrix by kmers
        features = svr_features_from_sequence(sequence, kmers)
        feature_size = len(features)
        if const_intercept: feature_size += 1 # If we are to use a const intercept term, we will have one more feature
        if model_dict['size'] != feature_size:
            raise Exception("Model size {} does not match feature size {}.\nPlease check paramaters for width, kmers, and const_intercept".format(model_size, feature_size))
        predictions, accuracy, values = predict(features, model_dict['model'], const_intercept)
        yield position, sequence, predictions[0]


if __name__ == '__main__':
    predict_genome('hg19.fa', 'GGAA', 36, 'ELK1_100nM_Bound_filtered_normalized_GGAA_1a2a3mer_format.model', [1,2,3], True, 'output.txt')
