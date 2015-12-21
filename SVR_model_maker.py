#!/usr/bin/env python

"""
Created on Nov 12, 2014

@author: josh
"""

import argparse
import itertools
import random
import string
import time
from distutils.util import strtobool
from operator import itemgetter
from subprocess import *

import numpy

# ===============================================================================
#  Command Line Arguments

parser = argparse.ArgumentParser(description='E2F Binding Model')
parser.add_argument('-i', metavar='PBMFile',
                    help='The results of a custom PBM experiment with sequences centered by binding site (i.e. using PWM)',
                    dest='pbmfile',
                    required=True)
parser.add_argument('-o', metavar='OutFilePrefix',
                    help='Optional, the prefix that all output files will be based on (do not include file extension). See program notes for proper format of this file',
                    dest='outprefix')
parser.add_argument('-g', "--gridsearch",
                    help="Flag for running a grid search, if optimal cost and epsilon values are not known",
                    action="store_true")
parser.add_argument('--seqlength', metavar='SequenceLength',
                    help="Change the length of the PBM sequence (Default is 36, new sequence will remain centered according to original 36mer from PBM data)",
                    dest="length")
parser.add_argument('--feature', metavar='FeatureType',
                    help="Define the type of features, i.e. 2 for 2mers, 123 for 1, 2, and 3-mers, etc; default is 3 for 3mers",
                    dest="rawkmer")
parser.add_argument('--extrafiles',
                    help="Print extra files: including all matrix files, feature definitions (sequence and position), and model sequence files",
                    action="store_true")
parser.add_argument('-c', metavar='SVR_cost',
                    help='The cost value input for LibSVM. If running a grid search, this should be a string of numbers in quotes, i.e. "0.05 0.1 0.5".',
                    dest='c')
parser.add_argument('-p', metavar='SVR_epsilon',
                    help='The epsilon value input for LibSVM. If running a grid search, this should be a string of numbers in quotes, i.e. "0.05 0.1 0.5".',
                    dest='p')
args = parser.parse_args()
pbmfile = args.pbmfile

### Getting the output prefix name
if args.outprefix:
    outprefix = args.outprefix
elif args.gridsearch:
    outprefix = os.path.splitext(pbmfile)[0] + '_SVR-gridsearch'
else:
    outprefix = os.path.splitext(pbmfile)[0] + '_SVR-model'

if args.extrafiles:
    extrafiles = 'yes'
else:
    extrafiles = 'no'

infofile = outprefix + '_run-info.txt'
f_info = open(infofile, 'a', 1)

# ===============================================================================

print "\n_____  Running SVR_model_maker _____"

''' Notes ======================================================================
PBM file format: Tab separated file has header line with column names, and
columns are: Name, ID, Sequence, Orientation1, Orientation2, Best-orientation,
and Replicate_intensity_difference.



'''

# pbmfile = 'E2F4-cust2-scores.txt'
# libsvm_generate_matrix_multimers_v5(pbmfile,36,3, 'good',12)
# testfile = 'E2F4-cust2-scores_regr-matrix-test_3mer-feat_36mer-seq_all_good-cores_0.5orientdiff_v5_linear_nocoresoutsidecenter12bp.txt'
# trainfile = 'E2F4-cust2-scores_regr-matrix-train_3mer-feat_36mer-seq_all_good-cores_0.5orientdiff_v5_linear_nocoresoutsidecenter12bp.txt'
# libsvm_run(0.1,0.15,trainfile,testfile)
# libsvm_feature_weights_V2(trainfile,3)


''' Called Parameters =========================================================='''
### These parameters are querried by the program

if args.gridsearch:
    if args.c:
        c_list_in = args.c
        c_list = c_list_in.split()
        c_list = [float(x) for x in c_list]
    else:
        print '\nEnter the different cost values to test, separated by spaces (i.e. "0.01 0.1 1")'
        c_list_in = raw_input(':')
        c_list = c_list_in.split()
        c_list = [float(x) for x in c_list]
    if args.p:
        p_list_in = args.p
        p_list = p_list_in.split()
        p_list = [float(x) for x in p_list]
    else:
        print '\nEnter the different epsilon values to test, separated by spaces'
        p_list_in = raw_input(':')
        p_list = p_list_in.split()
        p_list = [float(x) for x in p_list]
else:
    # Getting the libsvm cost and epsilon parameters
    if args.c:
        c = args.c
    else:
        print "\nLibSVM cost needs to be specified (if not known, try 0.01)"
        c = raw_input(':')
    if args.p:
        p = args.p
    else:
        print "\nLibSVM epsilon needs to be specified (if not known, try 0.01)"
        p = raw_input(':')

''' Other Parameters ================================================================='''
### Different parameters used during the process that can be changed from there defaults here

if args.length:
    length = int(args.length)
else:
    length = 36  # Assigning the default value
if args.rawkmer:
    rawkmer = args.rawkmer
else:
    rawkmer = "3"  # Assigning the default value

# Defining what kind of features we want
kmers = []
for item in str(rawkmer): kmers.append(int(item))
kinfo = ''
for x in kmers: kinfo = kinfo + str(x) + " + "
kinfo = kinfo[:-3] + " mer features"

# In the case of E2Fs, we want to make sure that good sequences have a good
# central core, without having a high-affinity 4-mer in the flanking sequences.
# Set this to [''] to use all sequences (i.e. every 4mer is good in the flanks
# and the core)
searchstrings = ['GCGG', 'CCGC', 'GCGC']

# Number bins we split the sequences into for SVR, where one bin is used for
# testing the model, and the remaining bins are combined for the training sequences.
svrbins = 5

''' Defining the modules ======================================================='''


def yes_no_query(question):
    sys.stdout.write('%s [y/n]\n' % question)
    while True:
        try:
            return strtobool(raw_input().lower())
        except ValueError:
            sys.stdout.write('Please respond with \'y\' or \'n\'.\n')


def read_data(filename):
    """ Creates an array for every cell in a tab separated text file"""
    data = []
    try:
        f = open(filename, 'r')  # opens the file as "f"
    except IOError:
        print "Could not open the file:", filename
        sys.exit()
    for line in f:  # for each line in the file
        l = string.split(line.strip(),
                         '\t')  # removes any carriage returns, then splits the tab separated line into columns
        data.append(l)  # Add each line to the array
    f.close()  # closing the file
    return data


def list_bins(l, bins):
    """Takes some list l and breaks it up into bins"""
    n = float(len(l)) / bins
    return [l[int(n * i):int(n * (i + 1))] for i in range(bins)]


def libsvm_generate_matrix(seqlist):
    """Generates the matrix file from a list of sequences and their scores"""
    svrmatrix = []
    for line in seqlist:
        score, seq = line
        ###Creating the list of features with relevent info
        featureinfo = [['feature', 'position', 'featnum', 'featvalue'],
                       ['start', 'na', 1, 1]]  # header, and first feature which never changes
        featnum = 2  # the first feature is already definde as 1:1, so we start with 2
        for k in kmers:
            ###Getting list of all possible combinates of bases, of length k, and all 4mers (for core)
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
                    featureinfo.append(
                            [feature, n1, featnum, featvalue])  # adding the info about the feature to the list
                    featnum += 1  # increasing the feature number by 1
        features = [score]  # starting a new list for building the matrix for SVR
        for x in range(1, len(featureinfo)):  # for every feature, in this list (skipping first item because header)
            features.append(
                    str(featureinfo[x][2]) + ':' + str(featureinfo[x][3]))  # putting the feature values into the list
        svrmatrix.append(
                features)  # adding the features for each sequence to the master list of features for this set of sequences
    return svrmatrix, featureinfo


def libsvm_run_gridsearch(p_list, c_list, pbmfile):
    """ Using libsvm, runs a grid search varying cost and epsilon, then prints a table of results, followed by the best R squared"""
    test = 0
    results = []

    ### Create your own module for selecting sequences from the PBM file for your protein
    print "\nFinding good sequences from the pbmfile to use for SVR"
    seqlist = read_pbm_sequences(pbmfile)
    check_data_length(seqlist)
    ### Creating the matrix file
    print "Generating matrix file for grid search..."
    svrmatrix = libsvm_generate_matrix(seqlist)[0]
    matrixfile = outprefix + '_gridsearch_matrix.txt'
    f = open(matrixfile, 'w')
    for line in svrmatrix: print >> f, ("\t".join(map(str, line)))
    f.close()

    ### Reading data in libsvm format
    print '\nDoing Libsvm grid search for', pbmfile, "with 5-fold cross validation"
    print >> f_info, '\nDoing Libsvm grid search for', pbmfile, "with 5-fold cross validation"

    for c in c_list:  # testing different values of C (cost)
        row = []
        pvals = [' ']
        for p in p_list:  # testing different values of epsilon, for each different value of C
            start_time = time.time()
            print "\nTesting epsilon =", p, "and cost =", c
            pvals.append(p)
            command = "svm-train -s 3 -v 5 -t 0 -c " + str(c) + " -p " + str(
                    p) + " " + matrixfile  # the command we want to run
            args = string.split(command)  # we need to split these into individual items
            output, error = Popen(args, stdout=PIPE,
                                  stderr=PIPE).communicate()  # running the command, and storing the results
            if len(error) > 0:  # if running the command caused an error, print the error, then quit the loop
                print error
                break
            out = string.split(output, '\n')  # Splitting the output by line
            SCC = float(string.split(out[-2])[
                            -1])  # the output is a list of strings (lines), the last word in the last line (with info) is the R squared
            print 'RSQ=', SCC
            print >> f_info, 'p=', p, ' c=', c, ' RSQ=', SCC
            print >> f_info, 'Round completed in', (time.time() - start_time) / 60, 'minutes'
            print 'Round completed in', (time.time() - start_time) / 60, 'minutes'
            if SCC > test:
                test = SCC
                best = [SCC, c, p]
            row.append(SCC)
        results.append([c] + row)
    results.insert(0, pvals)
    print >> f_info, '\nTable of results - columns are epsilon, and rows are cost\n'
    for line in results: print >> f_info, ("\t".join(map(str, line)))
    print >> f_info, '\nBest R squared is ', best[0], ' - with c = ', best[1], ' and p = ', best[2]
    print '\nBest R squared is ', best[0], ' - with c = ', best[1], ' and p = ', best[2]


PBM_SCORE_COLUMN = 5
PBM_SEQUENCE_COLUMN = 2


def read_pbm_sequences(pbmfile):
    data = read_data(pbmfile)  # Read the PBM file into a list of lists
    print "\nReading sequences from the PBM file"
    data = data[1:]  # Remove the header row
    # Extract the score and the sequence
    return [[float(row[PBM_SCORE_COLUMN]), row[PBM_SEQUENCE_COLUMN]] for row in data]


MAX_NUM_SEQUENCES = 5000


def check_data_length(seqlist):
    if len(seqlist) > MAX_NUM_SEQUENCES:
        print "Error: PBM file contains {} sequences, which exceeds the maximum of {}".format(len(seqlist),
                                                                                              MAX_NUM_SEQUENCES)
        sys.exit(1)


def libsvm_run(c, p, pbmfile):
    """ Using libsvm, for running the best set of values (best if obtained from a grid search), using the train and test matrix files"""

    ### Create your own module for selecting sequences from the PBM file for your protein
    print "\nFinding good sequences from the pbmfile to use for SVR"
    allseqlist = read_pbm_sequences(pbmfile)

    check_data_length(allseqlist)

    print "Using", len(allseqlist), "sequences"
    traincount = int(float(len(allseqlist)) * ((float(svrbins) - 1) / float(svrbins)))
    print "Using about", traincount, "sequences (+/-1) out of", len(
            allseqlist), "for training the model"  # plus or minus 1, for each of the binned sets
    print >> f_info, "Using about", traincount, "sequences (+/-1) out of", len(
            allseqlist), "for training the model"  # plus or minus 1, for each of the binned sets

    random.shuffle(allseqlist)  # randomizing the list
    seqbins = list_bins(allseqlist, svrbins)
    rsqinfo = []
    #### Takes the set of bins and uses one for the test sequences, and the rest for training the model
    for x in range(len(seqbins)):
        print "\nRunning SVR on round", x + 1, "of", len(seqbins)
        testseq = seqbins[x]
        trainseq = list(itertools.chain(*(
            seqbins[:x] + seqbins[
                          x + 1:])))  # Note, the remaining bins need to be collapsed into a single list of lists
        # print >>f_info, 'Number of sequences in the training data set is', len(trainseq), 'out of', len(allseqlist), 'sequences'
        print "Generating the matrix files..."
        trainmatrix, featureinfo = libsvm_generate_matrix(trainseq)
        testmatrix = libsvm_generate_matrix(testseq)[0]

        ### Writting matrix related output files
        trainmatrixfile = outprefix + '_train_matrix' + '_' + str(x + 1) + '.txt'
        f1 = open(trainmatrixfile, 'w', 1)
        for line in trainmatrix: print >> f1, ("\t".join(map(str, line)))
        f1.close()
        testmatrixfile = outprefix + '_test_matrix' + '_' + str(x + 1) + '.txt'
        f2 = open(testmatrixfile, 'w', 1)
        for line in testmatrix: print >> f2, ("\t".join(map(str, line)))
        f2.close()

        ###Training the model
        modelfile = trainmatrixfile + ".model"  # specify the model file name to svm-train, otherwise it writes in current dir
        print 'Training the model for run', x + 1, '...'
        command = "svm-train -s 3 -t 0 -c " + str(c) + " -p " + str(
                p) + " " + trainmatrixfile + " " + modelfile  # the command we want to run
        args = string.split(command)  # we need to split these into individual items
        output, error = Popen(args, stdout=PIPE,
                              stderr=PIPE).communicate()  # running the command, and storing the results
        if len(error) > 0:  # if running the command caused an error, print the error
            print error
        out = string.split(output, '\n')  # Splitting the output by line

        ###testing the model
        print 'Testing the model for run', x + 1, '...'
        outfile = testmatrixfile[0:-4] + '_SVR-prediction.txt'
        args = ["svm-predict", testmatrixfile, modelfile, outfile]
        output, error = Popen(args, stdout=PIPE,
                              stderr=PIPE).communicate()  # running the command, and storing the results
        if len(error) > 0:  # if running the command caused an error, print the error
            print error
        out = string.split(output, '\n')  # Splitting the output by line
        SCC = float(string.split(out[-2])[
                        -2])  # the output is a list of strings (lines), the last word in the last line (with info) is the R squared
        print 'Squared Correlation Coefficient for run', x + 1, 'is', SCC
        print >> f_info, 'Squared Correlation Coefficient for run', x + 1, 'is', SCC
        rsqinfo.append([SCC, x])

        ### Printing the extra files if we have the extrafiles flag
        if extrafiles == 'yes':
            if x == 0:  # We only need this file once, not for every iteration
                featurefile = outprefix + '_example-feature-list.txt'
                f0 = open(featurefile, 'w', 1)
                for line in featureinfo: print >> f0, ("\t".join(map(str, line)))
                f0.close()
        if SCC >= max([rsqinfo[i][0] for i in range(len(rsqinfo))]):
            ###organizing the results
            testdata = read_data(testmatrixfile)
            predictdata = read_data(outfile)
            if extrafiles == 'yes':
                trainseqfile = outprefix + '_train_sequences.txt'
                testseqfile = outprefix + '_test_sequences.txt'

    rsqinfo = sorted(rsqinfo, key=itemgetter(0), reverse=True)  # sorting the list by the first column

    bestrun = rsqinfo[0][1]
    rsqlist = [float(rsqinfo[i][0]) for i in
               range(len(rsqinfo))]  # getting just the r squared scores (first column) from the list
    mean_rsq = numpy.mean(rsqlist)  # getting the average r squared
    std_rsq = numpy.std(rsqlist)  # getting the standard deviation of the r squareds
    print "\nBest R squared is", rsqinfo[0][0], "from run", bestrun + 1
    print "Average R squared =", mean_rsq, "\nStandard deviation =", std_rsq,
    print >> f_info, "\nBest R squared is", rsqinfo[0][
        0], "from run", bestrun + 1, "- Mean = ", mean_rsq, "- Standard Deviation = ", std_rsq

    if extrafiles == 'yes':
        f3 = open(trainseqfile, 'w', 1)
        for line in trainseq: print >> f3, ("\t".join(map(str, line)))
        f3.close()
        f4 = open(testseqfile, 'w', 1)
        for line in testseq: print >> f4, ("\t".join(map(str, line)))
        f4.close()

    ### Keeping only the files for the best run, and renaming them
    for x in range(len(seqbins)):
        trainmatrixfile = outprefix + '_train_matrix' + '_' + str(x + 1) + '.txt'
        testmatrixfile = outprefix + '_test_matrix' + '_' + str(x + 1) + '.txt'
        modelfile = trainmatrixfile + ".model"
        outfile = testmatrixfile[0:-4] + '_SVR-prediction.txt'
        if x != bestrun:
            os.remove(trainmatrixfile)
            os.remove(testmatrixfile)
            os.remove(modelfile)
            os.remove(outfile)
        else:
            libsvm_feature_weights(modelfile)
            os.rename(trainmatrixfile, trainmatrixfile[:-6] + '.txt')
            os.rename(testmatrixfile, testmatrixfile[:-6] + '.txt')
            os.rename(modelfile, trainmatrixfile[:-6] + '.txt.model')
            os.rename(outfile, testmatrixfile[0:-6] + '_SVR-prediction.txt')

    bestresults = [['Actual-Intensity', 'Predicted-Intensity']]
    for n1 in range(len(testdata)):  # for every line, except the first, wich contains headders
        bestresults.append([testdata[n1][0], predictdata[n1][0]])
    resultsfile = outprefix + '_SVR-test_prediction-results.txt'
    f = open(resultsfile, 'w', 1)
    for line in bestresults: print >> f, ("\t".join(map(str, line)))
    f.close

    ### Writing info to the info file
    print >> f_info, "\nOutput files for best run are:"
    print >> f_info, ' ', outprefix + '_train_matrix.txt', "<-- The LibSVM matrix file for the sequences used to trian the model"
    print >> f_info, ' ', outprefix + '_test_matrix.txt', "<-- The LibSVM matrix file for the sequences used to test the model"
    if extrafiles == 'yes':
        print >> f_info, ' ', featurefile, "<-- The LibSVM matrix for the sequences used to trian the model"
        print >> f_info, ' ', outprefix + '_train_sequences.txt', "<-- The actual sequences and the corresponding scores used for training the model"
        print >> f_info, ' ', outprefix + '_test_sequences.txt', "<-- The actual sequences and the corresponding scores used for testing the model"
        print >> f_info, ' ', featurefile, "<-- List of definitions for the features (the feature sequence and it's position in the complete sequence"
    print >> f_info, ' ', outprefix + '_train_matrix.txt.model', '  <-- The model file that can be used by libsvm to predict binding affinities'
    print >> f_info, ' ', outprefix + '_test_matrix_SVR-prediction.txt', '  <-- The actual predicted intensity scores for the test set generated by LibSVM'
    print >> f_info, ' ', resultsfile, '  <-- Contains both the actual intensities, and predicted intensities for the test set'


def libsvm_feature_weights(modelfile):
    """ Getting feature weights, using the output from libsvm, with features of size k, works with sequences of variable sizes"""
    # print 'Getting feature weights from', modelfile
    model = read_data(modelfile)[6:]  # Starting with row 7, to avoid header lines
    k = kmers[-1]
    ### Creating a dictionary with the base/doublet for each feature
    bases = []
    for n in itertools.product('ACGT', repeat=k):
        bases.append(''.join(n))

    features = []
    featnum = 2  # starting the feature number at 2, because the first feature is just '1:1'
    #     print model[0][0]
    seqsize = (len(model[0][0].split()) - 2) / (len(bases)) + (
        k - 1)  # convoluted way of getting the size of the sequence used
    for x in range(1, seqsize + 1 - (k - 1)):  # for each position in the sequence, starting with 1 (instead of 0)
        for base in bases:
            # print featnum, "\t", x, "\t", base
            features.append([featnum, x, base])
            featnum += 1
            #         print features
    weights = {key: 0 for key in
               range(1, len(features) + 2)}  # dictionary to hold the value of the weights for each feature
    counts = {key: 0 for key in range(1, len(features) + 2)}
    #    print counts
    #    print len(weights), weights
    for n1 in range(len(model)):  # for every line in the svm model file
        line = string.split((model[n1][0]).strip())  # putting each element of the line into a list, called "line"
        #        print line
        for n2 in range(1, len(line)):  # for every feature (the first item is the sequence weight, not a feature)
            feature = string.split((line[n2]).strip(),
                                   ':')  # separating each element into the feature number, and the value
            weight = float(line[0]) * float(
                    feature[1])  # finding the weight of each feature (= sequence weight x feature value)
            weights[int(feature[0])] = weights[int(
                    feature[0])] + weight  # updating the dictionary to account for the total weights for each feature
            counts[int(feature[0])] += int(feature[1])
            #             print 'Feature number is', feature[0], 'and value is', feature[1], 'weight is', line[0], 'val is', weight, 'count is',counts[int(feature[0])]
            #    print weights
            #    for x in range(1,346): print (x+1), "\t", weights[x]
    output = [['Feature_Number', 'Feature_Weight', 'Position', 'Feature', 'Count']]
    print "Writing the output"
    for n3 in range(len(weights) - 1):
        #         print features[n3][0],weights[n3+2],features[n3][1],features[n3][2]
        output.append([features[n3][0], weights[n3 + 2], features[n3][1], features[n3][2], counts[n3 + 2]])
    # for line in output: print "\t".join(map(str,line))
    f = open(modelfile[:-10] + '-weights.txt', "w")
    for line in output: print >> f, "\t".join(map(str, line))
    f.close
    return output


''' Running the program ========================================================'''

### Testing to see if we need to normalize the data in the PBM file
pbmdata = read_data(pbmfile)
scores = [float(row[5]) for row in pbmdata[1:]]
maxval, minval = max(scores), min(scores)
if maxval > 1 or maxval < 0 or minval > 1 or minval < 0:  # if the scores are not between 0 and 1, we need to normalize everything
    print >> f_info, "Normalizing the data in the PBM file"
    for i in range(1, len(pbmdata)):
        for j in [3, 4, 5]:  # each of the columns we need to normalize
            score = float(pbmdata[i][j])
            pbmdata[i][j] = (score - minval) / (maxval - minval)
        pbmdata[i][6] = pbmdata[i][3] - pbmdata[i][4]  # getting the new difference
    newpbmfile = os.path.splitext(pbmfile)[0] + '_normalized.txt'
    f = open(newpbmfile, 'w')
    for line in pbmdata: print >> f, "\t".join(map(str, line))
    f.close()
    pbmfile = newpbmfile  # making sure we only use this new normalized file from now on

print >> f_info, '\n================================================================================\n', \
    '\nProgram started on', time.strftime("%m\%d\%Y %H:%M"), \
    '\nGenerating the SVR model file, using libsvm', \
    '\n\nParamaters used in this run:', \
    '\n ', " ".join(
        map(str, searchstrings)), '  <-- Good core 4-mers used for selecting sequences for building the model', \
    '\n ', args.pbmfile, '  <-- Input PBM file', \
    '\n ', outprefix, '  <-- Output file pfrefix', \
    '\n ', length, '  <-- Sequence Length', \
    '\n ', kinfo, '  <-- Feature type', \
    '\n ', svrbins, '  <-- Number bins we split the sequences into for SVR, where for each bin, it is used for testing, with the rest used for training the model', \
    '\n  Linear   <-- LibSVM support vector regression model type',
if args.gridsearch:
    print "Running a grid search. Use the results and re-run this program to refine another grid search, or run full LibSVM"
    print "\nCost values to be tested:", c_list, "\nEpsilon values to be tested:", p_list
    print >> f_info, '\n ', c_list, '  <-- LibSVM costs tested', '\n ', p_list, '  <-- LibSVM epsilons tested\n'
    libsvm_run_gridsearch(p_list, c_list, pbmfile)

else:
    print "Running full libsvm"
    print >> f_info, '\n ', c, '  <-- LibSVM cost', '\n ', p, '  <-- LibSVM epsilon\n'
    libsvm_run(c, p, pbmfile)

f_info.close()
