# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 14:13:10 2019

@author: HF
"""

import os
import wget

d = os.getcwd()
names = []
failed = []
#for i in os.listdir(d + '\\found_peps'):
#    names.append()
i = os.listdir(d+'\\found_peps')
lists = {}
current = 0
def get_url(pepname):
    print('Current: ' + pepname, end = '\r')
    global failed
    url = 'https://files.rcsb.org/download/' + pepname + '.pdb'
    try:
        wget.download(url, out = d + '\\pdb',  bar= None)
    except:
        newstring = pepname
        found = False
        while not found and len(newstring) != 0:
            newstring = newstring[:-1]
            print('                                                              ', end = '\r')
            print('Current: ' + newstring, end = '\r')
            url = 'https://files.rcsb.org/download/' + newstring + '.pdb'
            if (newstring + '.pdb') not in os.listdir(d + '\\pdb'):
                try:
                    wget.download(url, out = d + '\\pdb', bar= None)
                    found = True
                except:
                    found = False
        if not found:
            failed.append(pepname)

with open(d + '\\pep_titles.txt', 'w+') as f:
    f.write(i[0])
    get_url(i[0])
    for j in range(1,len(i)):
        f.write(',' + i[j])
        get_url(i[j])

with open(d + '\\failedpeps.txt', 'w+') as f:
    f.write(failed[0])
    for j in range(1,len(failed)):
        f.write(',' + failed[j])
