# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 15:25:27 2015

@author: timmonen
"""
import os

class Sample(object):
    
    _primer_dict = None
    _fragment_names = None
    _fragment_to_primer = None
    
    
    def __init__(self, samplename):
        self.name = samplename
        self.load_sample_info()
        self.get_table_filename()
        self.get_read_filenames()
        self.get_divide_summary_filenames()
        self.get_pre_map_filename()
        self.get_primer_positions()
        self.get_fragment_positions()
        self.get_fragment_trim_positions()
        self.get_fragment_output_names()
        self.get_data_foldername()
        self.get_trashed_read_names()

        
    @property
    def primer_dict(self):
        if self._primer_dict is None:
            from Bio.Seq import reverse_complement as rc
            self._primer_dict = {'F1': ['CTCAATAAAGCTTGCCTTGAGTGC', rc('ACTGTATCATCTGCTCCTGTRTCT')],
                   'F2': ['AAATTGCAGGGCYCCTAG', rc('CTRTTAGCTGCCCCATCTACATAG')],
                   'F3B': ['CACACTAATGATGTAARACARTTAACAG', rc('GGGATGTGTACTTCTGAACTTAYTYTTGG')],
                   'F4': ['CGGGTTTATTWCAGRGACAGCAGA', rc('GGGGTTAAYTTTACACATGGYTTTA')],
                   'F5a': ['GGCATYTCCTATGGCAGGAAGAAG', rc('GTGGTGCARATGAGTTTTCCAGAGCA')],
                   'F6': ['GGGTTCTTRGGARCAGCAGGAAG', rc('ATTGAGGCTTAAGCAGTGGGTTC')],}
        return self._primer_dict
        
    @property
    def fragment_names(self):
        if self._fragment_names is None:
            self._fragment_names = ["F1","F2","F3","F4","F5","F6"]
        return self._fragment_names
        
    @property
    def fragment_to_primer(self):
        if self._fragment_to_primer is None:
            self._fragment_to_primer = {'F1':'F1', 'F2':'F2','F3':'F3B','F4':'F4','F5':'F5a','F6':'F6',}
        return self._fragment_to_primer


    @staticmethod
    def write_json(data, filename):
        from .utils.formats import write_json as wj
        wj(data, filename)

    @staticmethod
    def get_data_foldername():
        import pipeline # is this just referring to the root?
        foldername = os.path.abspath(os.path.sep.join([pipeline.__path__[0],
                                                       os.path.pardir,
                                                       'data']))
        foldername = foldername+os.path.sep
        return foldername
        
        
    def get_primer_positions(self):     
        from Bio import SeqIO
        import seqanpy.seqanpy as sap
        import numpy as np
        ref_name = self.get_data_foldername()+'NL4-3.fasta'   
        ref = str(SeqIO.read(ref_name, "fasta").seq)  
        primer_positions = {}
        for n, prim in enumerate(self.primer_dict.keys()):
            # Take care of LTRs
            subref = ref[:7000] if prim[:2] in ['F1', 'F2'] else ref[2000:]
            offset = 0 if prim[:2] in ['F1', 'F2'] else 2000
            # Forward primer sequence
            prim_fwd =  self.primer_dict[prim][0]
            pos_fwd = np.array(list(sap.align_overlap(subref, prim_fwd)[2]))
            pos_fwd1 = np.argmax(pos_fwd != '-')
            pos_fwd2 = pos_fwd1 + len(prim_fwd)
            # FBackward primer sequence
            prim_bwd =  self.primer_dict[prim][1]
            pos_bwd = np.array(list(sap.align_overlap(subref, prim_bwd)[2]))
            pos_bwd1 = np.argmax(pos_bwd != '-')
            pos_bwd2 = pos_bwd1 + len(prim_bwd)
        
            primer_positions[prim] = [[pos_fwd1 + offset,pos_fwd2 + offset],[pos_bwd1 + offset,pos_bwd2 + offset]]
        return primer_positions
        
    def get_fragment_positions(self):
        fragment_positions = {}
        for n, frag in enumerate(self.fragment_names):
            fragment_positions[frag] = [self.get_primer_positions()[self.fragment_to_primer[frag]][0][0],
                               self.get_primer_positions()[self.fragment_to_primer[frag]][1][1]]
        return fragment_positions
        
        
    def get_fragment_trim_positions(self):
        fragment_trim_positions = {}
        for n, frag in enumerate(self.fragment_names):
            fragment_trim_positions[frag] = [self.get_primer_positions()[self.fragment_to_primer[frag]][0][1]+1,
                               self.get_primer_positions()[self.fragment_to_primer[frag]][1][0]-1]
        return fragment_trim_positions
                               
    def get_fragment_output_names(self):
        fragment_output_names = []
        for i, frag in enumerate(self.fragment_names):
            fragment_output_names.append(self.get_data_foldername() + frag + '.sam')
        return fragment_output_names
        
    def get_trashed_read_names(self):
        trashed_read_names = []
        trashed_reads = ['unmapped_read','small_insert', 'big_insert', 'trashed_primer', 'ambiguous_fragment', 'no_fragment']
        for i, trashed in enumerate(trashed_reads):
            trashed_read_names.append(self.get_data_foldername() + trashed + '.sam')
        return trashed_read_names
        
    def get_table_filename(self):
        foldername=self.get_data_foldername()
        filename = foldername+'sample_table.csv'
        return filename

        
    def load_sample_info(self):
        import pandas as pd
        filename = self.get_table_filename()
        table = pd.read_csv(filename, sep='\t', index_col=0)
        table.index = map(str, table.index)
        sample = table.loc[self.name]
        self.run = sample['seq run']
        self.adapter = sample['adapter']        

        
    def get_read_filenames(self, gzip=True, trimmed=False):
        '''Get the filenames of the demultiplexed reads'''
        foldername=self.get_data_foldername()
        #filenames = ['read1', 'read2']
        filenames = ['read1_10000', 'read2_10000']
        for i,fn in enumerate(filenames):
            fn = foldername+fn
            if trimmed:
                fn = fn+'_trimmed'
            fn = fn+'.fastq'
            if gzip:
                fn = fn+'.gz' 
            filenames[i] = fn
            
        if not trimmed:
            summary = None
        else:
            summary = foldername+'trim_summary.json'
            
        return {'data': filenames,
                'summary': summary}
                
    def get_divide_summary_filenames(self):
        '''Get the filenames of the demultiplexed reads'''
        foldername=self.get_data_foldername()
        summary = foldername+'trim_and_divide_summary.json'
        return summary
                
    
    def get_pre_map_filename(self):
        foldername=self.get_data_foldername()
        filename = foldername+'pre_map.sam'
        summary = foldername+'pre_map_summary.json'
        return {'data': filename,
                'summary': summary}
        
        
        
    def trim_reads(self, *args, **kwargs):
        from .trim_reads import trim_reads as tr
        return tr(self, *args, **kwargs)
        
                        
    def pre_map(self, *args, **kwargs):
        from .pre_map import pre_map as pm
        return pm(self, *args, **kwargs)
        

    def trim_and_divide(self, *args, **kwargs):
        from .trim_and_divide import trim_and_divide as tad
        return tad(self, *args, **kwargs)

