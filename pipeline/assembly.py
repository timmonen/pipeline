# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 11:40:13 2015

@author: timmonen"""

import sys
import os
import subprocess as sp
import argparse
import pysam
from pipeline.sample import Sample    

def reads_to_seqrecord(reads):
    '''Build a FASTQ record out of BAM reads
    
    Note: copied from Bio.SeqIO.QualityIO.py
    '''
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    # Precompute conversion table
    SANGER_SCORE_OFFSET = ord("!")
    q_mapping = dict()
    for letter in xrange(0, 255):
        q_mapping[chr(letter)] = letter - SANGER_SCORE_OFFSET
    
    seqs = []
    for read in reads:
        # Get the sequence first
        descr = read.qname
        id = read.qname
        name = id
        from Bio.Alphabet import IUPAC
        from Bio.Alphabet import IUPAC
        from Bio.Seq import reverse_complement as rc
        if not read.is_reverse:
            record = SeqRecord(Seq(read.seq, IUPAC.ambiguous_dna),
                           id=id, name=name, description=descr)
        if read.is_reverse:
            record = SeqRecord(Seq(rc(read.seq), IUPAC.ambiguous_dna),
                           id=id, name=name, description=descr)
    
        # Get the qualities second
        qualities = [q_mapping[letter] for letter in read.qual]
        if qualities and (min(qualities) < 0 or max(qualities) > 93):
            raise ValueError("Invalid character in quality string")
        dict.__setitem__(record._per_letter_annotations,
                         "phred_quality", qualities)

        seqs.append(record)

    return (seqs)
    

    


# Functions
def assembly(sample, VERBOSE=0,**kwargs):
    input_filenames = sample.get_fragment_output_names()
    data_folder = sample.get_data_foldername()
    from Bio import SeqIO
    for f in range(len(input_filenames)):
        with pysam.Samfile(input_filenames[f], 'r') as samfile:
            fastqreads = reads_to_seqrecord(samfile)
            SeqIO.write(fastqreads, data_folder + sample.fragment_names[f] + ".fastq",'fastq')
        sp.call(['/home/timmonen/.local/bin/iva', '--fr',data_folder + sample.fragment_names[f] + ".fastq",data_folder + sample.fragment_names[f]])

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
    sample.assembly(VERBOSE=args.verbose)

