#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 14:06:12 2019

@author: gemsec-user
"""

import pickle
import pdb_parser as pp
import os
d = os.getcwd()

# takes all of the cleaned pdb files and constructs peptides from them, dumps that into pickle into storage
def make_all_peps(length):
    maxV = pp.find_max(length)
    for pep in os.listdir(d + '/cleaned'):
        done = pp.make_structure(pep, length , maxV, cleaned = True) # done is an array of "finished" dictionaries from pdb parser
        # serializes all of the finished dictonaries to files
        for i in done:
            pickle.dump(i, open(d + '/made_adj/' + pep.replace('.pdb','') + '/' + i['sequence'] + '.pickle', 'wb'))   