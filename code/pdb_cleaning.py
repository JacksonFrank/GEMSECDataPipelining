# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 17:26:55 2019

@author: jtfl2

Helper script to clean pdb folders
"""

import cleaning
import os
import subprocess
import platform
import easygui

"""
This check is for my local machine, as graphical interfaces
    on Windows Subsystem for Linux is different
"""
if 'Microsoft' in platform.release():
    os.environ["DISPLAY"] = ":0"
    print "Display: "
    subprocess.call("echo $DISPLAY", shell=True)



d = os.getcwd()

easygui.fileopenbox('Choose a file', '', '')

for i in os.listdir('pdb/'):
    cleaning.cleanATOM(d + '/pdb/' + i, out_file= d + '/cleaned/' + i, ext = '.pdb')