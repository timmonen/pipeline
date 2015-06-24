"""
Spyder Editor

This temporary script file is located here:
/home/timmonen/.spyder2/.temp.py
"""
from Bio import SeqIO
from Bio import SeqRecord, Seq
import numpy as np

folder_name = "/home/timmonen/projects/example/"

reads1 = SeqIO.parse(folder_name+"read1.fastq", "fastq")

#Record object

reads2 = SeqIO.parse(folder_name+"read2.fastq", "fastq")

myfilename = folder_name+"read1_trim.fastq"
with open(myfilename, 'a') as f:
    for read in reads1:
        qold = np.array(read.letter_annotations['phred_quality'], int)
        ind = qold > 30
        seqold=np.array(read.seq)
        seqnew = seqold[ind].tostring()
        qnew = list(qold[ind])
        read.letter_annotations = {}
        read.seq=Seq.Seq(seqnew)
        read.letter_annotations['phred_quality']=qnew
        SeqIO.write(read,f,'fastq')
num=0
myfilename = folder_name+"read2_trim.fastq"
with open(myfilename, 'a') as f:
    for read in reads2:
        num=num+1
        print(num)
        qold = np.array(read.letter_annotations['phred_quality'], int)
        ind = qold > 30
        seqold=np.array(read.seq)
        seqnew = seqold[ind].tostring()
        qnew = list(qold[ind])
        read.letter_annotations = {}
        read.seq=Seq.Seq(seqnew)
        read.letter_annotations['phred_quality']=qnew
        SeqIO.write(read,f,'fastq')
        