# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 17:26:55 2019

@author: jtfl2

Helper script to clean pdb folders
"""

import cleaning
import os
import platform
import easygui

if 'Microsoft' in platform.release():
    os.system('export DISPLAY=:0')
    os.system('echo $DISPLAY')



d = os.getcwd()

easygui.fileopenbox('Choose a file', '', '')

for i in os.listdir('pdb/'):
    cleaning.cleanATOM(d + '/pdb/' + i, out_file= d + '/cleaned/' + i, ext = '.pdb')