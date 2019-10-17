#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 11:41:59 2019

@author: gemsec-user
"""
import numpy as np
import prody as pd
import PeptideBuilder as pb
import os
import Bio
import cleaning
from decimal import Decimal
#from pdbtools import pdbtools as pdb

ElementSymbols = 'H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe, Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn, Fr, Ra, Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr, Rf, Db, Sg, Bh, Hs, Mt, Uun, Uuu, Uub, Uut, Uuq, Uup, Uuh, Uus, Uuo'.split(', ')
aa = 'A,N,D,C,Q,E,G,H,I,L,K,M,F,P,R,S,T,W,Y,V'.split(',')
aa3 = 'ALA,ASN,ASP,CYS,GLN,GLU,GLY,HIS,ILE,LEU,LYS,MET,PHE,PRO,ARG,SER,THR,TRP,TYR,VAL'.split(',')
element_index = ['N', 'C', 'O', 'S', 'Se']
d = os.getcwd()
parsed_aa = {}
def parse_aa():
    global parsed_aa
    for amino in aa:
        if not os.path.exists(d + '/amino_acids'):
            os.mkdir(d + '/amino_acids')
        out = Bio.PDB.PDBIO()
        i = pb.make_structure(amino, [180]*len(amino),[180]*len(amino))
        out.set_structure(i)
        out.save(d + '/amino_acids/' + amino + '.pdb')
        cleaning.cleanATOM(d + '/amino_acids/' + amino + '.pdb', out_file= d + '/amino_acids/' + amino + '.pdb', ext = '.pdb')
        temp = pd.parsePDB(d + '/amino_acids/' + amino + ".pdb")
        parsed_aa[amino] = []
        for atom in temp.iterAtoms():
            parsed_aa[amino].append(str(atom.getName()))
        
def remove_padding(nodes):
    atoms = []
    current = 0
    col = nodes[0:5, current]
    while sum(col) != 0:
        atoms.append(element_index[col.tolist().index(1.0)])
        current += 1
        col = nodes[0:5, current]
    return atoms

    
def heuristic(index, node, amino_acid):
    correct = 0
    total = 0
    for atom in parsed_aa[amino_acid]:
        if (index+total) < len(node) and ElementSymbols[int(node[index+total][0]) - 1] == atom[0]:
            correct += 1
        total += 1
    return float(correct/total)
    
def find_sequence_recurs(nodes, length, current_ind, current_sequence, current_value):
    if len(parsed_aa.keys()) == 0:
        parse_aa()
    if len(current_sequence) == length:
        global POSSIBLE_SEQUENCES
        if current_value in POSSIBLE_SEQUENCES:
            POSSIBLE_SEQUENCES[current_value].append(current_sequence)
        else:
            POSSIBLE_SEQUENCES[current_value] = [current_sequence]
    values = []
    for a in aa:
        values.append(heuristic(current_ind,nodes, a))
    max_value = max(values)
    if max_value > 0.8:
        for i in range(len(values)):
            if max_value == values[i]:
                amino = aa[i]
                find_sequence_recurs(nodes, length, current_ind + len(parsed_aa[amino]), current_sequence + amino, current_value + max_value)
    
def find_white_space(total_space, text): 
    return ' '*(total_space - len(text))
        
POSSIBLE_SEQUENCES = None
def decode(encoding, save_loc = d, save_name = '', find_coord = False, use_coord = False):
    if len(parsed_aa.keys()) == 0:
        parse_aa()
    if save_name == '':
        save_name = encoding['sequence'] + '.pdb'
    placed = []
    new_nodes = remove_padding(encoding['index'])
    if not use_coord:
        D = encoding['secondary']
        placed.append([new_nodes[0], (0,0,0)])
        placed.append([new_nodes[1], (D[0,1],0,0)])
        x = (D[1,2]**2 - D[0,2]**2 - D[0,1]**2)/(-2 * D[0,1])
        y = (abs(D[0,2]**2 - x**2))**(0.5)
        placed.append([new_nodes[2], (x,y,0)])
        P = placed[2][1][0]**2 + placed[2][1][1]**2
        for i in range(3,len(new_nodes)):
            x = (D[1,i]**2 - D[0,i]**2 - D[0,1]**2)/(-2*D[0,1])
            y = (D[2,i]**2 - D[0,i]**2 - P + (2*x*placed[2][1][0]))/(-2*placed[2][1][1])
            z = (abs(D[0,i]**2 - x**2 - y**2))**(0.5)
            placed.append([new_nodes[i], (x,y,z)])
        if find_coord:
            final = np.zeros((len(encoding['secondary'][0]),3))
            for i in range(len(placed)):
                final[i, 0] = placed[i][1][0]
                final[i, 1] = placed[i][1][1]
                final[i, 2] = placed[i][1][2]
            return final
    else:
        for i in range(3,len(new_nodes)):
            placed.append([new_nodes[i], (encoding['coordinates'][i][0],encoding['coordinates'][i][1],encoding['coordinates'][i][2])])
    with open(save_loc + '/' + save_name, 'w+') as g:
        counter = 0
        amino_num = 0
        for i in range(len(placed)):
            if counter == 0:
                counter = len(parsed_aa[encoding['ele_to_amino'][i][1]])
                amino_num += 1
            string = 'ATOM' #+ str(i + 1) + '  '+ encoding['seq_to_atoms'][i][0]
            string += find_white_space(7, str(i + 1)) + str(i+1) + '  '
            string += encoding['ele_to_amino'][i][0] + find_white_space(4, encoding['ele_to_amino'][i][0])
            string += aa3[aa.index(encoding['ele_to_amino'][i][1])] + ' A'
            string += find_white_space(4, str(amino_num)) + str(amino_num)
            string += find_white_space(12, str(round(Decimal(placed[i][1][0]), 3))) + str(round(Decimal(placed[i][1][0]), 3))
            string += find_white_space(8, str(round(Decimal(placed[i][1][1]), 3))) + str(round(Decimal(placed[i][1][1]), 3))
            string += find_white_space(8, str(round(Decimal(placed[i][1][2]), 3))) + str(round(Decimal(placed[i][1][2]), 3))
            string += '  1.00  0.00'
            string += find_white_space(11, placed[i][0]) + placed[i][0]
            g.write(string + '\n')
            counter -= 1
    return save_loc + '/' + save_name
        
        