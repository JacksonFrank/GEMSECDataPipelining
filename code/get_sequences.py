#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 12:36:24 2019

@author: gemsec-user
"""

import os
import motif as m
import prody as pd
import pdb_parser as pp
#import pickle
bad = ['X', 'Z', 'B', 'O']

# For every pdb file in given directory, 
def sequence_maker(length, pdb_folder = os.getcwd() + '/pdb', outputname = 'allpeps.txt'):
    outputname = str(length) + outputname
    sequences = {}
    total = len(os.listdir(pdb_folder))
    count = 0
    # for every file in the pdb folder
    for i in os.listdir(pdb_folder):
        count += 1
        print(str(count) + ' out of '  + str(total) + ' completed, Unique Found: ' + str(len(sequences.keys())), end = '\r')
        # gets an array of Atom Objects
        file = pd.parsePDB(pdb_folder + '/' + i)
        if file is None:
            os.remove(pdb_folder + '/' + i)
            continue
        seq = pp.find_seq(file)[0]
        # gets the peps from the sequence
        info = m.pep_parser(seq, length)
        for mot in info:
            if mot != 'full':
                for pep in info[mot]:
                    if pep not in sequences:
                        sequences[pep] = 1
                    else:
                        sequences[pep] += 1
    print('Saving Sequences')
    with open(outputname, 'w+') as g:
        g.write('Sequence,Amount Found \n')
        for i in sequences:
            g.write(i + ',' + str(sequences[i]) + '\n')
            
#sequence_maker(7, pdb_folder = '/home/gemsec-user/Desktop/Jons_stuff/pytorch_stuff/cleaned')
            
# prints all of the unique sequences in the given file
def sequence_counter(filename):
    with open(filename, 'r') as g:
        unique = 0
        not_unique = 0
        string = g.readline()
        string = g.readline().replace('\n','').split(',')
        while len(string) > 1:
            unique += 1
            not_unique += int(string[1])
            print('Unique: ' + str(unique) + ' Not Unique: ' + str(not_unique), end= '\r')
            string = g.readline().replace('\n','').split(',')
        print()
            
sequence_counter('/home/gemsec-user/Desktop/Jons_stuff/pytorch_stuff/7allpeps.txt')
        
        
        
        
        
        
        
        