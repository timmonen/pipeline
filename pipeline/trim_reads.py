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
import re
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
    read_name = '@'+read[0]
    seq_orig= read[1]
    qual_store = read[2]
    
    qual_orig = np.fromstring(read[2], np.int8)-SANGER_SCORE_OFFSET
    qual_ind = qual_orig < quality
    qual_str= map(str,1*qual_ind)
    qual_str = ''.join(qual_str)
    
    #

    # Weed out very small errors in mid-sequence that we will ignore
    startinsert= [match.start() for match in re.finditer(re.escape("10"), qual_str)]
    endinsert= [match.start() for match in re.finditer(re.escape("01"), qual_str)]
    
    si=np.array(startinsert)
    ei=np.array(endinsert)
    
    # If every score perfect, just return as is!
    if len(si) == 0 and len(ei) == 0:
        seq_insert = seq_orig
        qual_insert = qual_store
        cut_off=[0,len(seq_insert)]  
        print "perfect"
    else:
        # If no beginning to first insert, make one at -1
        if len(si) > 0 and len(ei) > 0:
            if ei[0] < si[0]:
                si = np.insert(si,0,-1)
        if len(si) == 0 and len(ei) > 0:
            si = np.int_([-1])
        if len(si) > 0 and len(ei) == 0:
            if len(si) > 0 and len(ei) == 0:
                ei = np.int_([len(qual_str)-1])
        # If no end to last insert, make one at end of sequence
        if si[len(si)-1] > ei[len(ei)-1]:
            ei = np.insert(ei,len(ei),len(qual_str)-1)
        # Now we have lined up start-end insert pairs
        si_compare = si[1:]
        ei_compare = ei[:-1]
        # Can make it fancier here, don't drop if two 010 in a row (1010) or (010010)
        # Do this base
        drop_ind = (si_compare-ei_compare) == 1
        remove_singles = si_compare[drop_ind]
        
        qual_list = list(qual_str)
        for remove in iter(remove_singles):
            qual_list[remove]='0'
        qual_edit = ''.join(qual_list)
    
    # After filtering out sites with bad quality scores we choose to ignore, get insert
    
        startinsert= [match.start() for match in re.finditer(re.escape("10"), qual_edit)]
        endinsert= [match.start() for match in re.finditer(re.escape("01"), qual_edit)]
        si = np.array(startinsert)
        ei = np.array(endinsert)
        # If no beginning to first insert, make one at -1
        if len(si) == 0 and len(ei) == 0:
            seq_insert = seq_orig
            qual_insert = qual_store
            cut_off=[0,len(seq_insert)]  
            print "perfect"
        else: 
            if len(si) > 0 and len(ei) > 0:
                if ei[0] < si[0]:
                    si = np.insert(si,0,-1)
            if len(si) == 0 and len(ei) > 0:
                si = np.int_([-1])
            if len(si) > 0 and len(ei) == 0:
                ei = np.int_([len(qual_str)-1])
            # If no end to last insert, make one at end of sequence
            if si[len(si)-1] > ei[len(ei)-1]:
                ei = np.insert(ei,len(ei),len(qual_str)-1)
            # Length of all possible inserts
            all_inserts = ei-si
            max_insert = max(ei-si)
            # Choose biggest insert
            insert_ind = all_inserts == max_insert
            
            if len(insert_ind[insert_ind]) > 1:
                pick = np.random.random_integers(0,len(insert_ind[insert_ind])-1,1)
                seq_insert = seq_orig[si[insert_ind][pick]+1:ei[insert_ind][pick]+1]
                qual_insert = qual_store[si[insert_ind][pick]+1:ei[insert_ind][pick]+1]
                cut_off=[range(len(seq_orig)+1)[si[insert_ind][pick]+1],range(len(seq_orig)+1)[ei[insert_ind][pick]+1]]
                print("Sampled insert randomly from equally sized ones")
                
            else:
                seq_insert = seq_orig[si[insert_ind]+1:ei[insert_ind]+1]
                qual_insert = qual_store[si[insert_ind]+1:ei[insert_ind]+1]
                cut_off=[range(len(seq_orig)+1)[si[insert_ind]+1],range(len(seq_orig)+1)[ei[insert_ind]+1]]

    read_new = (read_name,seq_insert,qual_insert)

    return (read_new)

    
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
             #trims = [trim_read(read, **kwargs) for read in reads]
             trims= [trim_read(read) for read in reads]
             lrs = map(len,[trims[0][1],trims[1][1]])
             if (lrs[0] > minlen_read1) and (lrs[1] > minlen_read2):
             #Join list to string
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
    return n_good
    

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
    n_runs = sample.trim_reads(VERBOSE=args.verbose)
    
