# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 13:24:57 2024

@author: Omer
"""
import os, sys
reedSolomonProjectDir = os.environ.get('REEDSOLOMON')
if reedSolomonProjectDir == None: 
     reedSolomonProjectDir = "c:/users/omer/reedSolomon/reedSolomon/"
sys.path.insert(0, reedSolomonProjectDir)
projectDir = os.environ.get('IEEE8032DJ')
#You don't have to define an environment variable, but then you have to give a path to the project here
if projectDir == None: 
     projectDir = "c:/users/omer/802.3/"
sys.path.insert(1, projectDir)
from arithmetic import gf256
from ieeeConstants import parityMatrix_177_5, g_177_1, IEEE_8023_INT_DATA_TYPE
from ieee8023dj_d0p1 import encode_177_5
from hammingBasics import simpleHammingDecoder
import numpy as np




def check_equivalenctToPolynomial():
    
    alpha = gf256([0,0,0,0,0,0,1,0])
    # Wathc out !!! changing a class variable !!!
    #alpha.generatorPolynomial == polynomial([1,1,0,0,1,1,0,1,1])
    wan = gf256(1)
    temp = gf256([0,0,0,0,0,0,0,1])
    result = np.zeros((8,256), IEEE_8023_INT_DATA_TYPE)
    i = 0
    while temp != wan:
        result[:,i] = temp.coefficients
        i = i + 1
        temp = temp * alpha
    return result