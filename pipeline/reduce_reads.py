# vim: fdm=marker
'''
'''

# Modules
import sys
import gzip
from itertools import izip
import argparse
import numpy as np
from Bio import SeqIO

folder_name = "/home/timmonen/projects/example/data/"

reads1 = SeqIO.parse(folder_name+"read1.fastq", "fastq")

reads2 = SeqIO.parse(folder_name+"read2.fastq", "fastq")


num=0
myfilename = folder_name+"read1_100000.fastq.gz"
with gzip.open(myfilename, 'a') as f:
     for read in reads1:
        num=num+1
        print num

        SeqIO.write(read,f,'fastq')
        if num >= 100000:
            break
            
myfilename = folder_name+"read2_100000.fastq.gz"
with gzip.open(myfilename, 'a') as f:
    for num, read in enumerate(reads2):
        print num

        SeqIO.write(read,f,'fastq')
        if num + 1 >= 100000:
            break
        