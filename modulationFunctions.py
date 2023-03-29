# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 13:03:02 2023

@author: Megatron
"""
import os
import numpy as np
import time
from numba import jit, int32, float32, jitclass, types, typed, boolean, float64, int64
#import math
import wifiMatrices

projectDir = os.environ.get('8023')
if projectDir == None:
    import pathlib
    projectDir = pathlib.Path(__file__).parent.absolute()

import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, projectDir)
#import io

LDPC_LOCAL_PRNG = np.random.RandomState(7134066)
# LDPC_**_DATA_TYPE stores the data type over which all arithmetic is done.
# It is a nice way of changing the data type of the entire implementation at one place.
LDPC_DATA_TYPE = np.int64
LDPC_INT_DATA_TYPE = np.int64
LDPC_DECIMAL_DATA_TYPE = np.float64
LDPC_SEED_DATA_TYPE = np.int64
# Omer Sella: Major breakdown warning: the bool data type is used to create a mask. Replacing it with int32 breaks the decoder.
LDPC_BOOL_DATA_TYPE = boolean
# Omer Sella: seeds can be integers between 0 and 2**31 - 1
LDPC_MAX_SEED = 2**31 - 1

NUMBA_INT = int64
NUMBA_FLOAT = float64
NUMBA_BOOL = boolean



@jit(nopython = True)
def slicer(vector, length):
    ## Omer Sella: slicer puts a threshold, everything above 0 is translated to 1,  otherwise 0 (including equality). Do not confuse with the reserved function name slice !
    slicedVector = np.ones(length, dtype = LDPC_INT_DATA_TYPE)#LDPC_DECIMAL_DATA_TYPE)
    slicedVector[np.where(vector <= 0)] = 0
    return slicedVector


@jit(nopython = True)
def modulate(vector, length):
    modulatedVector = np.ones(length, dtype = LDPC_DECIMAL_DATA_TYPE)
    modulatedVector[np.where(vector == 0)] = -1
    return modulatedVector