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
currentPath = os.getcwd()
sys.path.insert(1, currentPath)
import numpy as np
from ieeeConstants import tv1_parity, tv1_tp4, generatorMatrix_177_5
from ieee8023dj_d0p1 import encode_177_5,  g_177_1
def test_encode_177_5_tv_1():
    tv1Encoded = encode_177_5(G = g_177_1, M = tv1_tp4)
    assert np.all(tv1Encoded[60:68] == tv1_parity)
    
def test_compareEncodingToTestVector():
    M_xor = np.array([ (tv1_tp4[2*k] + tv1_tp4[2 * k + 1]) %2 for k in range(60) ])
    parityAsCalculatedInDraft = M_xor.dot(g_177_1) %2
    assert np.all(parityAsCalculatedInDraft == tv1_parity)

def test_bchEncoder():
    pathToGeneratorMatrix = currentPath + "/bchMatrixEncoder.npy"
    h = np.load(pathToGeneratorMatrix)
    
#def test_checkEqual():
#    cw, mxor = encode_177_5(generatorMatrix_177_5, tv1_tp4)
#    cw2 = encodeUsingMatrixOnly(generatorMatrix_177_5, mxor[1,0:60].transpose())
#    return np.all(cw2 == mxor[1,:].transpose())