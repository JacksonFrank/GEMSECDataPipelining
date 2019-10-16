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
def make_all_peps(length):
    maxV = pp.find_max(length)
    for pep in os.listdir(d + '/cleaned'):
        done = pp.make_structure(pep, length , maxV, cleaned = True)
        for i in done:
            pickle.dump(i, open(d + '/made_adj/' + pep.replace('.pdb','') + '/' + i['sequence'] + '.pickle', 'wb'))   