# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 14:54:46 2019

@author: jtfl2
"""

import prody as pd
import os
import numpy as np
import cleaning
import motif as m
from functools import partial
import PeptideBuilder as pb
from symbols import *
#from Bio import Bio.PDB
#import Bio 

class PDBParser: 
    d = os.getcwd()

    def find_max(length):
        return 14*length


    def set_temp_value(val):
        for i in range(val):
            if not os.path.exists(d + '/temp_pdb/' + str(i)):
                os.mkdir(d + '/temp_pdb/' + str(i))

    def find_index(full_pep, pep):
        for i in range(len(full_pep)):
            if pep == full_pep[i:i+len(pep)]:
                return [i, i+len(pep)]
        return -1

    def seperate_pdb(loc, temp_value):
        temp = []
        current = 0
        for i in os.listdir(d + '/temp_pdb/'+str(temp_value)):
            os.remove(d + "/temp_pdb/" +str(temp_value) +"/" + i)
        new = True
        current_amino = 0
        skip_amino = -1
        with open(loc,'r') as g:
            string = g.readline()
            while string:
                if skip_amino != int(string[23:26]):
                    if current_amino != int(string[23:26]):
                        new = True
                        current_amino = int(string[23:26])
                    if string.startswith('ATOM') and string[-4] != 'H':
                        if new and string[-4] == 'N':
                            new = False
                            temp.append(string)
                        elif not new:
                            temp.append(string)
                        else:
                            skip_amino = current_amino
                    elif string.startswith('TER'):
                        with open(d + "/temp_pdb/" +str(temp_value) +'/temp' + str(current) + '.pdb', 'w+') as h:
                            for i in temp:
                                h.write(i)
                        current += 1
                        temp = []
                string = g.readline()
            

    def find_seq(file):
        sequence = {}
        sequence['full'] = ''
        current = -1
        found_elements = {}
        for i in range(len(file)):
            if file[i].getName() == 'N':
                current += 1
                sequence['full'] += file[i].getSequence() #AA[AA3.index(string[17:20])]
                sequence[current] = [i]
                found_elements[file[i].getElement()] = 1
            elif current in sequence:
                sequence[current].append(i)
                found_elements[file[i].getElement()] = 1
        return sequence, found_elements
                    
    #    sequences = []
    #    seq = ''
    #    global seq_to_pdb
    #    seq_to_pdb = {}
    #    temp = {}
    #    start_count = 1
    #    current = 0
    #    with open(loc, 'r') as g:
    #        string = g.readline()
    #        while string
    #            count = int(string[4:11])
    #            if string.startswith('ATOM'):
    #                if current != int(string[23:26]):
    #                    current = int(string[23:26])
    #                    seq += AA[AA3.index(string[17:20])]
    #                    temp[len(seq)] = [int(string[4:11])]
    #                else:
    #                    temp[len(seq)].append(int(string[4:11]))
    #            elif string.startswith('TER'):
    #                seq_to_pdb[string[21]] = temp.copy()
    #                temp = {}
    #                sequences.append([seq, (start_count, count)])
    #                seq = ''
    #                start_count = count + 1
    #            string = g.readline()
    #    return sequence


    def possible_bonds(i, pep, ind):
        possible = []
        possible.append(pep[i-ind[0]])
        if i != ind[0]:
            possible.append(pep[i -ind[0] -1])
        if i != ind[1] - 1:
            possible.append(pep[i - ind[0] +1])
        return possible
        
    class memoize(object):
        def __init__(self, func):
            self.func = func
            self.cache = {}

        def __call__(self, *args):
            if args in self.cache:
                return self.cache[args]
            else:
                result = self.func(*args)
                self.cache[args] = result
                return result

        def __get__(self, obj, objtype):
            return partial(self.__call__, obj)
        
        
    def clean_all(pdb_input, pdb_output):
        print(pdb_input, pdb_output)
        total = len(os.listdir(pdb_input))
        count = 0
        for i in os.listdir(pdb_input):
            cleaning.cleanATOM(pdb_input + '/' + i, out_file= pdb_output + '/' + i, ext = '.pdb')
            count += 1
            print(str(count) + ' out of ' + str(total) + ' completed', end = '\r')


    def make_structure(pdb_loc, length, maxValues, temp_value):
    #    pdb = pdb_loc.split('/')[-1]
        seperate_pdb(pdb_loc, temp_value)
        completed = []
        for i in os.listdir(d + "/temp_pdb/" +str(temp_value)):
            try:
                pdb_parsed = pd.parsePDB(d + "/temp_pdb/" +str(temp_value) +"/" + i)
            except:
                continue
            if pdb_parsed == None:
                continue
            seq, elements = find_seq(pdb_parsed)
            for i in elements:
                if i not in ELEMENT_INDEX:
                    return []
            info = m.pep_parser(seq['full'], length)
    #        if not os.path.exists(d + '/made_adj/' + pdb.replace('.pdb','')):
    #            os.mkdir(d + '/made_adj/' + pdb.replace('.pdb',''))
            for motif in info:
                if motif == 'full':
                    continue
                for pep in info[motif]:
                    completed.append(main_calc(pep, seq, maxValues, pdb_parsed))
        return completed

    def make_from_small_pdb(file_loc):
        pdb_parsed = pd.parsePDB(file_loc)
        seq, elements = find_seq(pdb_parsed)
        maxV = find_max(len(seq['full']))
        return main_calc([seq['full'], 0, len(seq)-1], seq, maxV, pdb_parsed)

    def main_calc(pep_info, seq_info, mV, file):
        Nterm = False
        pep = pep_info[0]
    #    print(pep)
        ind = [pep_info[1], pep_info[2]]
        new_pep = []
        for i in range(ind[0],ind[1]):
            for j in seq_info[i]:
                new_pep.append([file[j], i])

        lengths = np.zeros((mV, mV))
        index = np.zeros((len(ELEMENT_INDEX) + 1,mV))
        bonded = np.zeros((mV, mV))
        amino_acid_list = []
        for i in range(len(new_pep)):
            ele = float(ELEMENT_SYMBOLS.index(str(new_pep[i][0].getElement()))) + 1
            if ele == 6.0:
                index[5][i] = 0.0
            elif ele == 7.0 and Nterm:
                index[5][i] = 1.0
            elif ele == 8.0:
                index[5][i] = 2.0
            if Nterm == False and new_pep[i][0].getName() == 'N':
                Nterm = True
            index[ELEMENT_INDEX.index(new_pep[i][0].getElement())][i] = 1
            amino_acid_list.append([str(new_pep[i][0].getName()), str(new_pep[i][0].getSequence())])
            for j in range(len(new_pep)):
                if i != j and lengths[i][j] == 0:
                    lengths[i][j] = (float(pd.calcDistance(new_pep[i][0],new_pep[j][0])))
                    lengths[j][i] = lengths[i][j]
                if (lengths[i][j] < 1.9) and (new_pep[i][0] != new_pep[j][0]):
                    if new_pep[j][0].getSequence() in possible_bonds(new_pep[i][1], pep, ind):
                        bonded[i][j] = lengths[i][j]
                        bonded[j][i] = bonded[i][j]
        
        final = {}
        final['sequence'] = pep
        final['ele_to_amino'] = amino_acid_list
        final['index'] = index
        final['secondary'] = lengths
        final['primary'] = bonded
        return final

    #test = make_from_small_pdb('/home/jonathan/Desktop/MLPHHGA.pdb')
    #import pickle
    #done = make_structure(d + '/cleaned/1akj.pdb',7, find_max(7))
    #for i in done:
    #    current = 0
    #    while os.path.exists(d + '/made_adj/1akj/' + i['sequence'] + str(current) + '.pickle'):
    #        current += 1
    #    pickle.dump(i, open(d + '/made_adj/1akj/' + i['sequence'] + str(current) +  '.pickle', 'wb'))    
    #         
                
                
                
                
                
                
                
                
                
                
                