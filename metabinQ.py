### MetaBinQ ###

import os
import math
import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
sns.set()
import warnings
from Bio import SeqIO
from sklearn.decomposition import PCA
from sklearn.cluster import OPTICS, KMeans
from sklearn.mixture import GaussianMixture
# own python extension
import tnCounter


def check_minseqlen(value):
    ivalue = int(value)
    if ivalue < 2500:
      raise argparse.ArgumentTypeError("minseqlen < 2500 is not allowed")
    return ivalue

def check_windowsize(value):
    ivalue = int(value)
    if ivalue < 4000:
      raise argparse.ArgumentTypeError("windowsize < 4000 is not allowed")
    return ivalue

def options():
    parser = argparse.ArgumentParser()
    
    # user options
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-i", "--input", help = "Bins in fasta file format", nargs="+", default = [])
    group.add_argument("-c", "--counts", help = "Pre computed count files", nargs="+", default = [])
    parser.add_argument("-d", "--counts2disk", help = "Write count files to disk", action='store_true')
    parser.add_argument("-o", "--output", help = "Output file [default: metabinQ.out]", default="metabinQ.out")
    parser.add_argument("-f", "--force", help = "Force recalculation", action='store_true')
    parser.add_argument("-v", "--verbose", help = "Verbose output [default: false]", default = False, action='store_true')
    
    # developer options
    parser.add_argument("-m", "--minseqlen", help = argparse.SUPPRESS, default=2500, type=check_minseqlen)
    parser.add_argument("-w", "--windowsize", help = argparse.SUPPRESS, default=4000, type=check_windowsize)
    args = parser.parse_args()
    return(args)


def split_fragment(sequence, n):
    """
    Split sequence into fragments no longer than n
    
    Inputs:
    sequence - sequence
    n - max length of fragment
    
    Output - fragments
    """
    for i in range(0, len(sequence), n):
        yield sequence[i:i + n]
        

def split_overlap(sequence, n, overlap):
    """
    Split sequence into fragments no longer than n
    (with overlap)
    
    Inputs:
    sequence - sequence
    n - max length of fragment
    
    Output - fragments
    """
    for i in range(0, len(sequence), round(n/overlap)):
        yield str(sequence[i:i+n])


def opt_len(s,l):
    """add docstring"""
    return math.ceil(len(s)/round(len(s)/l))
    

def getContigTNF(path, min_length, windowsize, verbose=False):
    """add docstring"""
    contigs = SeqIO.parse(path, 'fasta')

    contigs_tnf = []
    for record in contigs:  
        if len(record.seq) > min_length:
            for fragment in split_overlap(record.seq, opt_len(record.seq, windowsize), 3):
                tnf = tnCounter.count(fragment)
                contigs_tnf.append(tnf)
        elif verbose:
            print("Warning: contig {} to small, it will be skipped".format(record.name)) 
    
    return(contigs_tnf)


def write_tnf(tnf, fname):
    """add docstring"""
    df = pd.DataFrame(tnf)
    df.to_csv(fname, header=True, index=False)


def new_z(contigs):
    mean = contigs.mean(axis=0)
    fractions = mean/mean.sum()
    E = np.dot(contigs.sum(axis=1).reshape(contigs.shape[0],1),
        fractions.reshape(1,contigs.shape[1]))
    z_scores = (contigs - E)/np.sqrt(E)
    return z_scores


def estimate_contamination(X, verbose):
    """ Add Docstring """

    pca = PCA()
    pca.fit(X)
    pca_data = pca.transform(X)

    clust = OPTICS(min_samples=50, xi=0.001, min_cluster_size=.05)
    
    if verbose:
        with warnings.catch_warnings():
          warnings.simplefilter("ignore")
          clust.fit(pca_data[:,:20])
    else:
        clust.fit(pca_data[:,:20])

    clusters = set(clust.labels_)
    if -1 in clusters:
        clusters.remove(-1)
    number_of_clusters = len(clusters)

    gmm = GaussianMixture(n_components=number_of_clusters, covariance_type='full')

    tmp = []
    for c in clusters:
        tmp.append(pca_data[clust.labels_ == c].mean(axis=0))

    est_contamination = 0 
    if number_of_clusters > 1:
        initial = np.array(tmp)

        gmm.fit(initial[:,:8])
        gmm.set_params(warm_start=True)

        labels = gmm.fit_predict(pca_data[:,:8])
        _, cluster_counts = np.unique(labels, return_counts=True)

        cont = np.min(cluster_counts)
        est_contamination = cont/np.sum(cluster_counts)

    return est_contamination

if __name__ == '__main__':

    args = options()

    if args.input and args.counts2disk:
      tmpdir = "tmp"
      if os.path.isdir(tmpdir):
          if args.force:
              print("Warning: tmp directory already exists and will be overwritten")
          else: 
            print("Error: count2disk option was set but tmp directory already exists".format(args.output))
            exit(0)
      else:
        os.mkdir(tmpdir)
   
    if os.path.isfile(args.output):
        if args.force:
            print("Warning: outputfile '{}' already exists and will be overwritten".format(args.output))
        else:
            print("Nothing to do, outputfile '{}' already exists".format(args.output))
            exit(0)

    outputfile = open(args.output, 'w')

    if args.input: # start from bins

      for f in args.input:
        if not os.path.isfile(f):
            print("Warning: Could not find file {}, it will be skipped".format(f))
            continue
        basename = os.path.basename(os.path.splitext(f)[0])
        
        # step1: calculate tetranculeotdie frequencies
        tnf = getContigTNF(f, args.minseqlen, args.windowsize, args.verbose)
        
        if args.counts2disk: # write count files if wanted
            fname = tmpdir + '/' + basename + ".csv"
            write_tnf(tnf, fname)
        
        # step2: analysis
        X = new_z(np.array(tnf))
        if len(tnf) > 150:
          contamination = estimate_contamination(X, args.verbose)
          outputfile.write("{},{}\n".format(basename, contamination))

    else: # start from count files

      for f in args.counts:
        if not os.path.isfile(f):
            print("Warning: Could not find file {}, it will be skipped".format(f))
            continue
        basename = os.path.basename(os.path.splitext(f)[0])

        # step1: read tetranculeotdie frequencies
        data = pd.read_csv(f)
        tnf = data.values

        # step2: analysis
        X = new_z(tnf)
        if len(tnf) > 150:
          contamination = estimate_contamination(X, args.verbose)
          outputfile.write("{},{}\n".format(basename, contamination))

    outputfile.close()

