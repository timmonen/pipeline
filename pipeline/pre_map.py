# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 17:48:33 2015

@author: timmonen
"""
import sys
import os
import subprocess as sp
import argparse
import pysam
from pipeline.sample import Sample    


# Functions
def pre_map_reads(sample):
    fn_in = sample.get_read_filenames(gzip=True, trimmed=True)['data']
    fn_outd = sample.get_pre_map_filename()
    fn_out = fn_outd['data']
    fn_outs = fn_outd['summary']
    
    sp.call(['/home/timmonen/bin/yara_indexer', sample.get_data_foldername()+'NL4-3.fasta', '-o', sample.get_data_foldername()+'NL4-3.index'])
    '''sp.call(['/home/timmonen/bin/yara_mapper', '-ll', '500', '-le', '300',
             '-o', fn_out,
             sample.get_data_foldername()+'NL4-3.fasta',
             fn_in[0], fn_in[1],
            ])'''
    sp.call(['/home/timmonen/bin/yara_mapper',
             sample.get_data_foldername()+'NL4-3.index',
             fn_in[0], fn_in[1], '-o', fn_out, '--ll', '500', '--le', '300'])

    seqs = []            
    with pysam.Samfile(fn_out, 'r') as bamfile:
        for ir, read in enumerate(bamfile):
            seqs.append(read.seq)
            if ir > 5:
                break
            pos_read = 0
            pos_ref = read.pos
            for (bt, bl) in read.cigar:
                if bt == 1:
                    pos_read += bl
                    # we ignore insertions
                    continue
                
                if bt == 2:
                    pos_ref += bl
                    # we ignore deletions
                    continue
                
                if bt == 0:
                    pos_read += bl
                    pos_ref += bl
                    # do something
                
                else:
                    continue
                




  # Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Pre map reads')
    parser.add_argument('--sample', required=True,
                        help='MiSeq sample to analyze')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()

    sample = Sample(args.sample)
    
    sample.pre_map_reads()
    