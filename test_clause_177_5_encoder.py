# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 15:30:10 2024

@author: Omer
"""
import os
import sys
projectDir = os.environ.get('IEEE8032DJ')
#You don't have to define an environment variable, but then you have to give a path to the project here
if projectDir == None: 
     projectDir = "c:/users/omer/802.3/"
sys.path.insert(1, projectDir)
import numpy as np
import ieeeConstants
from ieee8023dj_d0p1 import encode_177_5
def test_encode_177_5_tv_1():
    tv1Encoded, tv1xored = encode_177_5(G = ieeeConstants.G, M = ieeeConstants.tv1_doubled)
    assert np.all(tv1Encoded[0,120:128] == ieeeConstants.tv1_parity)
    

def test_bchEncoder():
    pathToGeneratorMatrix = projectDir + "/bchMatrixEncoder.npy"
    h = np.load(pathToGeneratorMatrix)
    