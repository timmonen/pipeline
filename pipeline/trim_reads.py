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
import argparse
import numpy as np
from Bio import SeqIO


# Functions
def trim_reads(samplename, VERBOSE=0, quality=25, blocksize=10,
               minlen_read1=100, minlen_read2=50):
    '''Trim low quality at the end of reads'''
    from pipeline.sample import Sample    
    
    sample = Sample(samplename)
    print sample.name
    print sample.run, sample.adapter
    
    fn_in = sample.get_read_filenames(gzip=True)
    fn_out = sample.get_read_filenames(gzip=True, trimmed=True)

    sys.exit()

    n_good = 0
    n_discarded = 0

    with gzip.open(fn_in[0], 'rb') as fin1, \
         gzip.open(fn_in[1], 'rb') as fin2, \
         gzip.open(fn_out[0], 'wb') as fout1, \
         gzip.open(fn_out[1], 'wb') as fout2:

        it1 = SeqIO.read(fin1, 'fastq')
        it2 = SeqIO.read(fin2, 'fastq')
        for irp, reads in enumerate(izip(it1, it2)):

            if VERBOSE >= 2:
                if not ((irp + 1) % 10000):
                    print irp + 1

            # Trim both reads
            trims = [trim_read(read, quality=quality, blocksize=blocksize)
                     for read in reads]

            lrs = map(len, trims)
            if (lrs[0] > minlen_read1) and (lrs[1] > minlen_read2):
                SeqIO.write(trims[0], fout1, 'fastq')
                SeqIO.write(trims[1], fout2, 'fastq')
                n_good += 1
            else:
                n_discarded += 1

    if VERBOSE:
        print 'Trim lowq ends of reads:'
        print 'Good:', n_good
        print 'Discarded:', n_discarded



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Trim loq quality end of reads')
    parser.add_argument('--sample', required=True,
                        help='Sample to analyze')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    # args contains sample and verbose argument
    args = parser.parse_args()
    trim_reads(args.sample, VERBOSE=args.verbose)
