# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 17:26:55 2019

@author: jtfl2

Cleans pdb folders
"""

import cleaning
import os
import subprocess
import platform
import easygui

def cleanPDB():
    """
    This check is for my local machine, as graphical interfaces
        on Windows Subsystem for Linux are different
    """
    if 'Microsoft' in platform.release():
        os.environ["DISPLAY"] = ":0"
        print "Display: "
        subprocess.call("echo $DISPLAY", shell=True)

    d = os.getcwd()

    filepath = easygui.fileopenbox('Choose a file', filetypes=["*.pdb"])

    if filepath is not None:
        cleaning.cleanATOM(filepath)