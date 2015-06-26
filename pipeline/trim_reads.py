# vim: fdm=marker
'''
author:     Fabio Zanini
date:       29/06/14
content:    This script is an optional step before premapping, to trim the low-q
            end of reads in order to improve premapping (i.e. we recover more
            reads).
'''
# Modules
import sys
import gzip
from itertools import izip
import numpy as np
from Bio import SeqIO
from Bio import Seq
import argparse

from pipeline.sample import Sample    

# Functions

def trim_read(read, quality=30, blocksize=3):
    '''Trim low quality at the end of a given read when at least blocksize 
       low quality scores in a row'''
     # Converts illumina character to PHRED SCORE
    SANGER_SCORE_OFFSET = ord("!")
    qual_orig = np.fromstring(read[2], np.int8)-SANGER_SCORE_OFFSET
    ind = qual_orig < quality
    s = map(str,1*ind)
    s = ''.join(s)
    substr = '1'*blocksize
    cut=s.find(substr)
    if cut == -1:
        pass
    else:
        seq_orig=np.array(read.seq)
        seq_trim = seq_orig[:cut].tostring()
        qual_trim = list(qual_orig[:cut])
        read.letter_annotations = {}
        read.seq=Seq.Seq(seq_trim)
        read.letter_annotations['phred_quality']=qual_trim
    return read
    
def trim_reads(sample, VERBOSE=0, minlen_read1=100, minlen_read2=50,
               **kwargs):
    '''Trim low quality at the end of reads to improve initial rough mapping'''
    from Bio.SeqIO.QualityIO import FastqGeneralIterator as FGI

    fn_in = sample.get_read_filenames(gzip=True)['data']
    fn_outdict=sample.get_read_filenames(gzip=True, trimmed=True)
    fn_outd = fn_outdict['data']
    fn_outs = fn_outdict['summary']

    n_good = 0
    n_discarded = 0

    with gzip.open(fn_in[0], 'rb') as fin1, \
         gzip.open(fn_in[1], 'rb') as fin2, \
         gzip.open(fn_outd[0], 'wb') as fout1, \
         gzip.open(fn_outd[1], 'wb') as fout2:
         for irp, reads in enumerate(izip(FGI(fin1), FGI(fin2))):
            # Trim both reads
             trims = [trim_read(read, **kwargs) for read in reads]
             lrs = map(len, trims)
             if (lrs[0] > minlen_read1) and (lrs[1] > minlen_read2):
             # Join list to string
             # This makes the fastq format of 4 lines, \n means new line
                 fout1.write('\n'.join([trims[0][0], trims[0][1], '+', trims[0][2]])+'\n')
                 fout2.write('\n'.join([trims[1][0], trims[1][1], '+', trims[1][2]])+'\n')      
                 n_good += 1
             else:
                 n_discarded += 1
    if VERBOSE:
        print 'Trim lowq ends of reads:'
        print 'Good:', n_good
        print 'Discarded:', n_discarded
        
    # Write summary file
    summary = {'sample name': sample.name,
               'sequencing run': sample.run,
               'adapter': sample.adapter,
               'number of read pairs': irp + 1,
               'number of discarded': n_discarded,
               'number of good': n_good,
               'histotgram of number of trimmed bases': np.zeros(301).tolist(),
              }    
    sample.write_json(summary, fn_outs)
             


    # USING SeqIO.read gives following error:
# raise ValueError("More than one record found in handle")

  # Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Trim low quality end of reads')
    parser.add_argument('--sample', required=True,
                        help='MiSeq sample to analyze')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()

    sample = Sample(args.sample)
    sample.trim_reads(VERBOSE=args.verbose)
    
