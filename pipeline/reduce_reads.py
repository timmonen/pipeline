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

reads11 = SeqIO.parse(folder_name+"read1_1000.fastq.gz", "fastq")

reads2 = SeqIO.parse(folder_name+"read2.fastq", "fastq")

myfilename = folder_name+"read1_1000.fastq.gz"

num=0
myfilename = folder_name+"read1_1000.fastq.gz"
with gzip.open(myfilename, 'a') as f:
     for read in reads1:
        num=num+1
        print num
        SeqIO.write(read,f,'fastq')
        if num >= 1000:
            break
            
myfilename = folder_name+"read2_1000.fastq.gz"
with gzip.open(myfilename, 'a') as f:
    for num, read in enumerate(reads2):
        print num
        SeqIO.write(read,f,'fastq')
        if num + 1 >= 1000:
            break

        