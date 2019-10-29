# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 16:13:15 2019

@author: jtfl2
"""

import pdb_parser as pp
import pickle
import os
import shutil
import numpy
import easygui as eg

d = os.getcwd()
current = 1
m = pp.find_max(7)
lists = [list(x) for x in numpy.array_split(numpy.array(os.listdir(d + '/cleaned')),12)]
current_list = int(eg.enterbox(msg = 'Which Temp?'))
failedlist = []

def make_struc(j, prot, current):
    done = pp.make_structure(d + '/cleaned/' + j,7, m, current_list)
    if len(done) > 0:
        if not os.path.exists('/media/jonathan/drive' + str(current) + '/peptides/' + prot):
            os.mkdir('/media/jonathan/drive' + str(current) + '/peptides/' + prot)
        for i in done:
            current_n = 0
            while os.path.exists('/media/jonathan/drive' + str(current) + '/peptides/' + prot + '/' + i['sequence'] + str(current_n) + '.pickle'):
                current_n += 1
            pickle.dump(i, open('/media/jonathan/drive' + str(current) + '/peptides/' + prot + '/' + i['sequence'] + str(current_n) +  '.pickle', 'wb'))
    else:
        failedlist.append(j)

if __name__ == '__main__':
    for ind, j in enumerate(lists[current_list]):
        print('Protein: ' + j, 'Finished: ' + str(ind) + '/' +str(len(lists[current_list])), 'Failed: ' + str(len(failedlist)), 'Current Chunk: ' + str(current_list), 'Current Drive: ' + str(current), end = '\r')
        prot = j.replace('.pdb','')
        stop = False
        for check in range(3):
            if os.path.exists('/media/jonathan/drive' + str(check + 1) + '/peptides/' + prot):
                stop = True
                break
        if not stop:
            total, used, free = shutil.disk_usage('/media/jonathan/drive' + str(current))
            if (used/total) > 0.93:
                current += 1
#            try:
            make_struc(str(j),prot,current)
#            except:
#                failedlist.append(j)
    print('Protein: ' + j, 'Finished: ' + str(ind) + '/' +str(len(lists[current_list])), 'Failed: ' + str(len(failedlist)), 'Current Chunk: ' + str(current_list))
    with open('failed' + str(current_list)+ '.txt','w+') as g:
        g.write(str(failedlist))