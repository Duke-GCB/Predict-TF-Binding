#!/usr/bin/env python

import argparse
import itertools
from svmutil import *
from math import exp
from multiprocessing import Process

from Bio import SeqIO, Seq

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
                print("Warning: sub-sequence '{}' not found in features".format(sub_seq))
            svr_features.extend(exploded)
    return svr_features


def generate_matching_sequences(sequence, core, width, core_start=None):
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

    # If sequence and core are strings, build BioPython sequences out of them
    # so they can be reverse-complemented
    if isinstance(core, str):
        core = Seq.Seq(core)
    if isinstance(sequence, str):
        sequence = Seq.Seq(sequence)
    sequence_width = len(sequence)
    core_width = len(core)
    # Need to search for core and reverse complement in the window region
    # If RC is found in the region, return the reverse-complement of the window instead
    # Also, if core is palindromic, need to return both regions and return best score
    core_rc = core.reverse_complement()
    max_start = sequence_width - width
    # The core positions are calculated relative to the window (and not the overall sequence)
    # This works as long as both core and width are same parity
    if core_start is None:
        core_start = (width - core_width) / 2
    for start in range(max_start + 1):
        end = start + width
        window_sequence = sequence[start:end]
        # If any of the bases in the window are unknown, we cannot predict on the sequence
        if 'N' in window_sequence:
            continue
        window_core = window_sequence[core_start:core_start + core_width]
        # If core is palindromic, return two sequences and let the caller decide which to use
        if core == core_rc and window_core == core:
            yield start, (str(window_sequence), str(window_sequence.reverse_complement()),)
        elif window_core == core:
            yield start, (window_sequence,)
        elif window_core == core_rc:
            yield start, (str(window_sequence.reverse_complement()),)

def read_fasta_idx(fasta_file):
    """
    Reads a fasta-formatted sequence file and returns a Bio.SeqIO.index
    :param fasta_file: File containing the sequences, such as hg19.fa
    :return: A Bio.SeqIO index
    """
    return SeqIO.index(fasta_file, 'fasta')


def get_sequence_named(idx, name):
    """
    Gets the sequence for a name
    :param idx: The BioPython sequence index
    :param name: sequence name
    :return: The sequence, uppercased
    """
    record = idx[name]
    return record.seq.upper()


def load_model(model_file, check_size=True):
    """
    Loads a svm model from a file and computes its size
    :param model_file: The file name of the model to load
    :return: A dictionary with keys model, file, and size
    """
    model = svm_load_model(model_file)
    model_dict = {'model': model}
    if check_size:
        model_dict['size'] = len(model.get_SV()[0]) - 1 # sv includes a -1 term that is not present in the model file, so subtract 1
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
    predictions = svm_predict([1.0], [svm_matrix], model, '-q')
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
    print(chrom, position, position + width, score, file=file_handle)


def predict_fasta(fasta_file, sequence_names, core, width, model_file, kmers, const_intercept,
                  check_size, core_start, transform_scores, output_file):
    """
    Generate predictions on the provided fasta file.
    Predictions will only be generated on sequences of width 'width', matching the nucleotides of 'core' in the center

    :param fasta_file: File name containing sequences, user-generated or a whole genome
    :param core: sequence of nucleotides to find when generating predictions
    :param sequence_names: List of sequence namnes in the fasta on which to predict (e.g. chroms)
    :param width: width, in bases, of the window on which to generate predictions
    :param model_file: Name of the svm model file to load
    :param kmers: List of integers (e.g. [1,2,3]) for base combination in prediction generation. Must match model generation parameters
    :param const_intercept: true or false - whether or not to add a constant term to the matrix generation. Must match model generation
    :param check_size: true or false - if true, check the model's size against the expected size given kmers/const_intercept/width
    :param core_start: Integer or None. If not None, used instead of width-len(core) / 2 for core search in window.
    :param transform_scores: true or false - whether or not to transform the scores using the transform_score function
    :param output_file: Output file to write, in bed format
    :return: None
    """

    print('Loading model', model_file)
    model_dict = load_model(model_file, check_size)

    print('Loading fasta', fasta_file)
    idx = read_fasta_idx(fasta_file)

    # If no sequences are named, use all in the index
    if sequence_names is None:
        sequence_names = sorted(idx.keys())

    # Iterate over sequence
    with open(output_file, 'w') as output:
        for sequence_name in sequence_names:
            print('Predicting on', sequence_name)
            for position, sequence, score in predict_sequence(idx, sequence_name, core, width, model_dict, kmers, const_intercept, core_start, transform_scores):
                print_bed(output, sequence_name, position, width, score)
    print('Done')


def predict_sequence(sequence_idx, sequence_name, core, width, model_dict, kmers, const_intercept, core_start, transform_scores):
    """
    Generates predictions for a single sequence
    :param sequence_idx: The indexed SeqIO object
    :param sequence_name: The name of a sequence (must be the id in the index)
    :param core: sequence of nucleotides to find when generating predictions
    :param width: width, in bases, of the window on which to generate predictions
    :param model_dict: Dictionary containing a loaded svm model at 'model' and optionally it's size at 'size'
    :param kmers: List of integers (e.g. [1,2,3]) for base combination in prediction generation. Must match model generation parameters
    :param const_intercept: true or false - whether or not to add a constant term to the matrix generation. Must match model generation
    :param core_start: Integer or None. If not none, used instead of width-len(core) / 2 for core search in window.
    :param transform_scores: true or false - whether or not to transform the scores using the transform_score function
    :return: A generator, yielding the start position, sequence, and prediction score
    """

    sequence = get_sequence_named(sequence_idx, sequence_name)
    print("Generating matching sequences for core {}, width {}".format(core, width))

    for position, matching_sequences in generate_matching_sequences(sequence, core, width, core_start):
        # generator returns a position, and a tuple of 1 or 2 sequences
        # If two sequences are returned, core is palindromic and can bind on either strand
        # So generate predictions for both and return the best
        # 4. Translate the sequences into SVR matrix by kmers
        best_prediction = None
        best_match = None
        for matching_sequence in matching_sequences:
            features = svr_features_from_sequence(matching_sequence, kmers)
            feature_size = len(features)
            if const_intercept: feature_size += 1 # If we are to use a const intercept term, we will have one more feature
            if 'size' in model_dict and model_dict['size'] != feature_size:
                raise Exception('Model size {} does not match feature size {}.\nPlease check parameters for width, '
                                'kmers, and const_intercept'.format(model_dict['size'], feature_size))
            predictions, accuracy, values = predict(features, model_dict['model'], const_intercept)
            if best_prediction is None or predictions[0] > best_prediction:
                best_prediction = predictions[0]
                best_match = matching_sequence
        if best_prediction is None:
            continue
        if transform_scores:
            best_prediction = transform_score(best_prediction)
        yield position, best_match, best_prediction


def predictable_chroms():
    """
    Returns the list of chromosomes on which we can predict: chr1, chr2...chr22, chrX, chrY, chrM
    :return: List of chromosome names
    """
    chroms = list()
    for i in range(1,23):
        chroms.append('chr' + str(i))
    for i in ('X','Y'):
        chroms.append('chr' + i)
    return chroms


def transform_score(score):
    # f(x) = 1 / ( 1 + exp(-x) )  to obtain only values between 0 and 1.
    return 1.0 / (1.0 + exp(0.0 - score))

def main():
    parser = argparse.ArgumentParser(description = 'TF Predictions Generator')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-g', metavar='GenomeFile',
                        help='Genome File in fasta format, containing sequence for each chrom as a separate entry',
                        dest='genome_fasta_file')
    group.add_argument('-s', metavar='SequenceFile',
                       help='Sequence in fasta format on which to predict',
                       dest='sequence_fasta_file')
    parser.add_argument('--chroms', '--seqnames', metavar='Chroms',
                        help='Optional List of chromosomes/sequence names from the sequence file to predict',
                        dest='sequence_names',
                        nargs='*',)
    parser.add_argument('-m', metavar='ModelFile',
                        help='The .model file generated from LibSVM, matching the specified core and kmers' ,
                        dest='model_files',
                        required=True,
                        nargs='+')
    parser.add_argument('-c', metavar='Core',
                        help='Sequence of nucleotides to search as center of predictions',
                        dest='cores',
                        required=True,
                        nargs='+')
    parser.add_argument('-w', metavar='Width',
                        type=int,
                        help='Width, in bases, of the window on which to generate predictions',
                        dest='width',
                        required=True
                        )
    parser.add_argument('-k', metavar='Kmers',
                        help='List of integers (e.g. 1,2,3) for combination of bases in prediction. Must match model',
                        dest='kmers',
                        type=int,
                        nargs='+',
                        required=True)
    parser.add_argument('-i',
                        action='store_true',
                        help='Whether or not to include the constant term in matrix generation. Must match model',
                        dest='const_intercept')
    parser.add_argument('--skip-size-check',
                        action='store_true',
                        help='Skip checking model size against kmers/width/const_interceipt. Size checking ensures model and parameters are recommended, but takes a few seconds on startup.',
                        dest='skip_size_check',
                        default=False)
    parser.add_argument('--core-start', metavar='Offset',
                        help='Override position to look for the core inside the width. Defaults to (width - core_width) / 2 if not specififed.',
                        type=int,
                        dest='core_start')
    parser.add_argument('-t',
                        action='store_true',
                        help='Transforms predictions with logistic function f(x) = 1 / ( 1 + exp(-x) )',
                        dest='transform_scores')
    parser.add_argument('-o', metavar='OutputFile',
                        help='Output file to write, in bed format',
                        dest='output_files',
                        nargs='*',
                        required=False)
    args = parser.parse_args()
    const_intercept = args.const_intercept or False
    check_size = not args.skip_size_check

    # Sequence and Genome are mutually exclusive
    # But selecting a genome should restrict to subset of predictable chroms
    if args.genome_fasta_file:
        fasta_file = args.genome_fasta_file
        sequence_names = args.sequence_names or predictable_chroms()
    else:
        fasta_file = args.sequence_fasta_file
        sequence_names = args.sequence_names

    output_files = args.output_files

    if not output_files:
        # no output files specified, determine from model/core
        output_files = ['predict_output_{}_core_{}.bed'.format(i + 1, core) for i, core in enumerate(args.cores)]

    # Ensure that model/core/output files are the same length
    if not len(output_files) == len(args.cores) == len(args.model_files):
        raise Exception('When specifying multiple cores/models/output files, you must specify the same number of each')

    if len(output_files) != len(set(output_files)):
        raise Exception('When specifying multiple output files, they must be unique')

    processes = list()
    for (model_file, core, output_file) in zip(args.model_files, args.cores, output_files):
        p = Process(target=predict_fasta, args=(fasta_file, sequence_names, core, args.width, model_file, args.kmers,
                                                const_intercept, check_size, args.core_start, args.transform_scores,
                                                output_file))
        processes.append(p)
        p.start()
    for p in processes:
        p.join()


if __name__ == '__main__':
    main()
