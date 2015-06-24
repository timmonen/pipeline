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


    @staticmethod
    def get_table_filename():
        import pipeline
        foldername = os.path.abspath(os.path.sep.join([pipeline.__path__[0],
                                                       os.path.pardir,
                                                       'data']))
        filename = foldername+os.path.sep+'sample_table.csv'
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
        return 'ciao'