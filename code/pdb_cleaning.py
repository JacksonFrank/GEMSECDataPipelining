# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 17:26:55 2019

@author: jtfl2, Jackson Frank

Cleans pdb folders
"""

import cleaning
import os
import easygui

# cleans the pdb file to work with Rosetta
class PDBCleaning:
    def cleanPDB():
        d = os.getcwd()
        filepath = easygui.fileopenbox('Choose a file', filetypes=["*.pdb"])
        if filepath is not None:
            cleaning.cleanATOM(filepath)