# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 17:28:24 2015

@author: timmonen
"""

def write_json(data, filename):
    import json
    with open(filename, 'w') as f:
        json.dump(data, f, indent=4)