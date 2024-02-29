# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 13:03:02 2023

@author: Megatron
"""
import os
import numpy as np
import time
#from numba import jit, int32, float32, types, typed, boolean, float64, int64
#import math
projectDir = os.environ.get('8023')
if projectDir == None:
    import pathlib
    projectDir = pathlib.Path(__file__).parent.absolute()
import sys
from ieeeConstants import *



# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, projectDir)
#import io

# LDPC_**_DATA_TYPE stores the data type over which all arithmetic is done.
# It is a nice way of changing the data type of the entire implementation at one place.


#@jit(nopython = True)
def slicer(vector):
    ## Omer Sella: slicer puts a threshold, everything above 0 is translated to 1,  otherwise 0 (including equality). Do not confuse with the reserved function name slice !
    length = vector.shape[0]
    slicedVector = np.ones(length, dtype = IEEE_8023_INT_DATA_TYPE)
    slicedVector[np.where(vector <= 0)] = 0
    return slicedVector


#@jit(nopython = True)
def modulatePAM2(vector):
    modulatedVector = np.ones(vector.shape[0], dtype = IEEE_8023_DECIMAL_DATA_TYPE)
    modulatedVector[np.where(vector == 0)] = -1
    return modulatedVector

def modulatePAM4(vector, greyCoded = True):
    #vector is assumed to be a binary vector of even length
    
    # note that the following line means the stream goes [MSB0,LSB0,MSB1,LSB1 ...]
    newVector = 2 * vector[0::2]  +  vector[1::2]
    modulatedVector = np.zeros(vector.shape[0] // 2, dtype = IEEE_8023_DECIMAL_DATA_TYPE)
    if not greyCoded:
        #00 --> -1
        #01 --> -1/3
        #10 --> 1/3 
        #11 --> 1
        
        modulatedVector[np.where(newVector == 0)] = PAM4_LEVEL_LOW
        modulatedVector[np.where(newVector == 1)] = PAM4_LEVEL_MID_LOW
        modulatedVector[np.where(newVector == 2)] = PAM4_LEVEL_MID_HIGH
        modulatedVector[np.where(newVector == 3)] = PAM4_LEVEL_HIGH
    else:
        #00 --> -1
        #01 --> -1/3
        #10 --> 1 
        #11 --> 1/3
         modulatedVector[np.where(newVector == 0)] = PAM4_LEVEL_LOW
         modulatedVector[np.where(newVector == 1)] = PAM4_LEVEL_MID_LOW
         modulatedVector[np.where(newVector == 2)] = PAM4_LEVEL_HIGH
         modulatedVector[np.where(newVector == 3)] = PAM4_LEVEL_MID_HIGH
    return modulatedVector

def pam4Slicer(vector, greyCoded = True):
    pam4Symbols = np.zeros((vector.shape[0]), dtype = IEEE_8023_INT_DATA_TYPE)
    bitsDemodulated = np.zeros((2 * vector.shape[0]), dtype = IEEE_8023_INT_DATA_TYPE)
    if greyCoded == True:
        pam4Symbols[np.where(vector <= -(2/3))]                         = 0
        mask = (vector > -(2/3)) & (vector <=  0  )
        pam4Symbols[np.where(mask)]                                     = 1
        mask = (vector > 0) & (vector <=  (2/3))
        pam4Symbols[np.where(mask)]                                     = 3
        pam4Symbols[np.where( (vector > (2/3))) ]                       = 2
    else:
        pam4Symbols[np.where(vector <= -(2/3))]                         = 0
        mask = (vector > -(2/3)) & (vector <= 0)
        pam4Symbols[np.where(mask)]                                     = 1
        mask = (vector > 0) & (vector <=  2/3  )
        pam4Symbols[np.where(mask)]                                     = 2
        pam4Symbols[np.where( (vector > 2/3)) ]                         = 3
    for i in range(vector.shape[0]):
        symbol = pam4Symbols[i]
        #if greyCoded:
            #print(symbol)
        if symbol == 0:
            bitsDemodulated[2 * i ]    = 0
            bitsDemodulated[2 * i + 1] = 0
        elif symbol == 1:
            bitsDemodulated[2 * i ]    = 0
            bitsDemodulated[2 * i + 1] = 1
        elif symbol == 2:
            bitsDemodulated[2 * i ]    = 1
            bitsDemodulated[2 * i + 1] = 0
        else:
            bitsDemodulated[2 * i ]    = 1
            bitsDemodulated[2 * i + 1] = 1
    return bitsDemodulated, pam4Symbols