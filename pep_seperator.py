# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 13:58:14 2019

@author: jtfl2
"""

import os
import easygui as eg
import motif as m
import pickle

d = os.getcwd()
file = eg.fileopenbox(msg = 'Where is the file?')
peps = {}
full_string = ''
with open(file, 'r') as g:
    string = g.read()
    while string:
        full_string += string
        string = g.read()
full_string.replace('  ',' ')
full_string = full_string.split('>')
for i in range(len(full_string)):
    full_string[i] = full_string[i].replace('\n',' ').replace('>','').split(' ')
    if len(full_string[i]) > 1:
        peps[full_string[i][0]] = full_string[i][-2]
  
count = 0
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
        