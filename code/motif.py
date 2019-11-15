# -*- coding: utf-8 -*-
"""
Created on Tue May 21 16:22:23 2019

@author: jtfl2
"""
import os
import pickle
from symbols import *
d = os.getcwd()
functions = ["ϕX{2,3}ϕX{2,3}ϕXϕ","YXXϕ","P[TS]AP","QXXϕXX[FHT][FHY]","[ST]X[IL]P",
             "[RK]XLX{0,1}ϕ","PXIX[IV]","RXXLXXϕ","[ST]PXX[ST]","ϕ(K)X[DE]","[DE]X([ST])ϕ",
             "[DE]XXD[GSAN]","[ST](P)"]

bulky_hydrophobic = 'W Y V L F I'.split(' ')

functions = [x.replace('p','').replace(',','').replace('\t','').replace('|','') for x in functions]
combos = {}
#
#def make_combo(l, current):
#    new_list = []
#    if len(current) != 0:
#        for k in current:
#            for a in l:
#                new_list.append(k + a)
#    else:
#        for a in l:
#            new_list.append(a)
#    return new_list
#
#def possible_combinations():
#    global combos
#    for i in functions:
#        if i == "ϕX{23}ϕX{23}ϕXϕ":
#            combos[i] = [9,11]
#            continue
#        print('Working on ' + i)
#        c = []
#        minv = 1000
#        maxv = 0
#        cont = 0
#        for ind,j in enumerate(i):
#            print('current: [' + j + ',' + str(ind) + ']', end = '\r')
#            if cont != 0:
#                cont -= 1
#                continue
#            if j == 'X':
#                c = make_combo(aa, c)
#            elif j == 'ϕ':
#                c = make_combo(bulky_hydrophobic, c)
#            elif j == '[':
#                new_j = [ind + 1,i[ind  + 1]]
#                aa_list = []
#                while new_j[1] != ']':
#                    aa_list.append(new_j[1])
#                    new_j[0] += 1
#                    new_j[1] = i[new_j[0]]
#                cont = len(aa_list)
#                c = make_combo(aa_list, c)
#            elif j == '{':
#                num_list = []
#                new_j = [ind + 1, i[ind + 1]]
#                while new_j[1] != '}':
#                    num_list.append(int(new_j[1]))
#                    new_j[0] += 1
#                    new_j[1] = i[new_j[0]]
#                copied = False
#                c_copy = c.copy()
#                for k in num_list:
#                    values = c_copy.copy()
#                    for l in range(k):
#                        values = make_combo(aa, values)
#                    if not copied:
#                        c = values
#                        copied = True
#                    else:
#                        for v in values:
#                            c.append(v)
#                        
#            elif j in aa:
#                if len(c) != 0 :
#                    for k in range(len(c)):
#                        c[k] += j
#                else:
#                    c.append(j)
#        if not os.path.isdir(d + '\\motifs\\' + i):
#            os.mkdir(d + '\\motifs\\' + i)
#        with open(d + '\\motifs\\' + i + '\\combos.txt', 'w+') as f:
#            first = True
#            for comb in c:
#                minv = min(minv, len(comb))
#                maxv = max(maxv, len(comb))
#                if first:
#                    f.write(comb)
#                    first = False
#                f.write(',' + comb)
#        combos[i] = [minv,maxv]
#    with open('combo.pickle', 'wb') as f:
#        pickle.dump(combos, f, pickle.HIGHEST_PROTOCOL)
                    
motifs = {}                
    

# need more info on arguments
def make_peps(location, pep, length, pos = False):
    if length > (location[1] -  location[0]): 
        index = location[1] - length
        smaller = False
    else:
        index = location[0] - (length-1)
        smaller = True
    if index < 0:
        index = 0
    found = []
    while (not smaller and index != location[0]) or (smaller and index != location[1]):
        if (index + length) < len(pep):
            new = pep[index: index + length]
            found.append((new, index, index + length))
        index += 1
    return found

# need more info on arguments
def fun_tester(t, pep):
    for seq in t:
        found = 0
        for i in range(len(pep)):
            if seq[i] == 'X':
                if pep[i] in aa:
                    found += 1
            elif seq[i] == 'ϕ':
                if pep[i] in bulky_hydrophobic:
                    found += 1
            else:
                if seq[i] == pep[i]:
                    found += 1
        if found == len(seq):
            return True
    return False

# what calls these functions?
def fun1(pep):
    t = ['ϕXXϕXXϕXϕ']
    if len(pep) == 9:
        t = ['ϕXXϕXXϕXϕ']
    elif len(pep) == 10:
        t = ['ϕXXXϕXXϕXϕ','ϕXXϕXXXϕXϕ']
    elif len(pep) == 11:
        t = ['ϕXXXϕXXXϕXϕ']
    return fun_tester(t,pep)
motifs[functions[0]] = fun1
combos[functions[0]] = [9,11]

def fun2(pep):
    t = ['YXXϕ']
    return fun_tester(t,pep)
motifs[functions[1]] = fun2
combos[functions[1]] = [4,4]

def fun3(pep):
    t = ['PTAP','PSAP']
    return fun_tester(t,pep)
motifs[functions[2]] = fun3
combos[functions[2]] = [4,4]

def fun4(pep):
    t = ["QXXϕXXFF", "QXXϕXXFH", "QXXϕXXFY", "QXXϕXXHF", "QXXϕXXHH", "QXXϕXXHY", "QXXϕXXYF", "QXXϕXXYH", "QXXϕXXYY"]
    return fun_tester(t,pep)
motifs[functions[3]] = fun4
combos[functions[3]] = [8,8]

def fun5(pep):
    t = ["SXIP", "SXLP","TXIP","TXLP"]
    return fun_tester(t,pep)
motifs[functions[4]] = fun5
combos[functions[4]] = [4,4]

def fun6(pep):
    if len(pep) == 4:
        t = ["RXLϕ","KXLϕ"]
    else:
        t = ["RXLXϕ","KXLXϕ"]
    return fun_tester(t,pep)
motifs[functions[5]] = fun6
combos[functions[5]] = [4,5]

def fun7(pep):
    t = ["PXIXI","PXIXV"]
    return fun_tester(t,pep)
motifs[functions[6]] = fun7
combos[functions[6]] = [5,5]

def fun8(pep):
    t = ['RXXLXXϕ']
    return fun_tester(t,pep)
motifs[functions[7]] = fun8
combos[functions[7]] = [7,7]

def fun9(pep):
    t = ['SPXXS','TPXXS','SPXXT','TPXXT']
    return fun_tester(t,pep)
motifs[functions[8]] = fun9
combos[functions[8]] = [5,5]

def fun10(pep):
    t = ['ϕKXD', 'ϕKXE']
    return fun_tester(t,pep)
motifs[functions[9]] = fun10
combos[functions[9]] = [4,4]

def fun11(pep):
    t = ["DXSϕ", "DXTϕ", "EXSϕ", "EXTϕ"]
    return fun_tester(t,pep)
motifs[functions[10]] = fun11
combos[functions[10]] = [4,4]

def fun12(pep):
    t = ["DXXDG", "DXXDS", "DXXDA", "DXXDN", "EXXDG", "EXXDS", "EXXDA", "EXXDN"]
    return fun_tester(t,pep)
motifs[functions[11]] = fun12
combos[functions[11]] = [5,5]

def fun13(pep):
    t = ["SP","TP"]
    return fun_tester(t,pep)
motifs[functions[12]] = fun13
combos[functions[12]] = [2,2]

# where are combos supposed to be defined?
def pep_parser(pep, length, pos = False):
    found_motifs = {}
    found_motifs['full'] = pep
    for i in combos:
#        print()
        found_motifs[i] = []
        for j in range(len(pep)):
#            print('Current Motif: ' + i + '     Position: ' + str(j) + '/' + str(len(pep)) + '    Found = ' + str(len(found_motifs[i])), end = '\r')
            for size in range(combos[i][0], combos[i][1] + 1):
                if motifs[i](pep[j:j+size]):
                    l = make_peps((j,j+size),pep,length, pos)
                    for x in l:
                        found_motifs[i].append(x)
    return found_motifs


    


        
        
#if __name__ == "__main__":       
#    if os.path.exists(d + '/combo.pickle'):
#        with open('combo.pickle', 'rb') as f:
#            # The protocol version used is detected automatically, so we do not
#            # have to specify it.
#            combos = pickle.load(f)
#    else:
#        possible_combinations()
#peps = pep_parser('YSPRTPDRVSETDIQRLLHGVMEQLGIARPRVEYPAHQAMNLVGPQSIEGGAHEGLQHLGPFGNIPNIVAELTGDNTPKDFSEDQGYPDPPNPCPIGKTDDGCLENTPDTAEFSREFQLHQHLFDPEHDYPGLGKWNKKLLYEKMKGGQRRKRRSVNPYLQGQRLDNVVAKKSVPHFSDEDKDPE',7, pos= True)
#    print(peps)
