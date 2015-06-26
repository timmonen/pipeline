"""
Spyder Editor

This temporary script file is located here:
/home/timmonen/.spyder2/.temp.py
"""
from Bio import SeqIO
from Bio import SeqRecord, Seq
import numpy as np

folder_name = "/home/timmonen/projects/example/data/"

reads1 = SeqIO.parse(folder_name+"read1_1000.fastq", "fastq")

myfilename = folder_name+"read1_trim.fastq"
with open(myfilename, 'a') as f:
    for read in reads1:
        qual_orig = np.array(read.letter_annotations['phred_quality'], int)
        ind = qual_orig < 30
        s = map(str,1*ind)
        s = ''.join(s)
        substring = '1'*error_length
        cut=s.index(substring)
        seq_orig=np.array(read.seq)
        seq_trim = seq_orig[:cut].tostring()
        qual_trim = list(qual_orig[:cut])
        read.letter_annotations = {}
        read.seq=Seq.Seq(seq_trim)
        read.letter_annotations['phred_quality']=qual_trim
        SeqIO.write(read,f,'fastq')

# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Trim low quality end of reads')
    parser.add_argument('--sample', required=True,
                        help='MiSeq sample to analyze')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
