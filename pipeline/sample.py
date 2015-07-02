# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 15:25:27 2015

@author: timmonen
"""
import os

class Sample(object):
    def __init__(self, samplename):
        self.name = samplename
        self.load_sample_info()
        self.get_table_filename()
        self.get_read_filenames()
        
    @property
    def primer_dict(self):
        '''Give a dictionary of the primer sequences for each fragment'''
        from Bio.Seq import reverse_complement as rc
        primers = {'F1': ['CTCAATAAAGCTTGCCTTGAGTGC', rc('ACTGTATCATCTGCTCCTGTRTCT')],
                   'F2': ['AAATTGCAGGGCYCCTAG', rc('CTRTTAGCTGCCCCATCTACATAG')],
                   'F3B': ['CACACTAATGATGTAARACARTTAACAG', rc('GGGATGTGTACTTCTGAACTTAYTYTTGG')],
                   'F4': ['CGGGTTTATTWCAGRGACAGCAGA', rc('GGGGTTAAYTTTACACATGGYTTTA')],
                   'F5a': ['GGCATYTCCTATGGCAGGAAGAAG', rc('GTGGTGCARATGAGTTTTCCAGAGCA')],
                   'F6': ['GGGTTCTTRGGARCAGCAGGAAG', rc('ATTGAGGCTTAAGCAGTGGGTTC')],}
        return primers
        
    
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
        filenames = ['read1_1000', 'read2_1000']
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
                
    
    def get_pre_map_filename(self):
        foldername=self.get_data_foldername()
        filename = foldername+'pre_map.sam'
        summary = foldername+'pre_map_summary.json'
        return {'data': filename,
                'summary': summary}
        
        
        
    def trim_reads(self, *args, **kwargs):
        from .trim_reads import trim_reads as tr
        return tr(self, *args, **kwargs)
        
                        
    def pre_map_reads(self, *args, **kwargs):
        from .pre_map import pre_map_reads as pmr
        return pmr(self, *args, **kwargs)


