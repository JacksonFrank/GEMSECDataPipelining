# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 13:58:14 2019

@author: jtfl2
"""

import os
import easygui as eg
import motif as m
import pickle

# d is current directory
d = os.getcwd()
# o
file = eg.fileopenbox(msg = 'Where is the file?')
# peps are peptides
peps = {}
full_string = ''
# gets entire pdb seqres file as single string
with open(file, 'r') as g:
    string = g.read()
    while string:
        full_string += string
        string = g.read()
# replaces all double spaces with single spaces
full_string.replace('  ',' ')
# splits string into an array of proteins
full_string = full_string.split('>')
for i in range(len(full_string)):
    # splits the string into its components
    full_string[i] = full_string[i].replace('\n',' ').replace('>','').split(' ')
    if len(full_string[i]) > 1:
        # index of -2 means second to last element in the array
        # keys each (peptide? the first part of the seqres file) to the name of the peptide (I'm not sure this will work correctly, look at file formatting)
        peps[full_string[i][0]] = full_string[i][-2]
count = 0
# writes peptide information to file system
for i in peps:
    count += 1
    with open('pep_titles.txt','w+') as f:
        f.write(i + ' ')
    print('Current Peptide: ' + i + " (" +str(count) +"/"+str(len(peps)) + ")")
    if not os.path.isdir(d + '\\found_peps\\' + i):
        os.mkdir(d + '\\found_peps\\' + i)
    found_peps = m.pep_parser(peps[i],12)
    with open(d + '\\found_peps\\' + i + '\\peps.pickle', 'wb') as f:
        pickle.dump(found_peps, f, pickle.HIGHEST_PROTOCOL)
        