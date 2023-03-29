# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 13:05:28 2023

@author: Megatron
"""
import os
import numpy as np
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

@jit()
def addAWGN(vector, length, SNRdb, prng):
    ## The input SNR is in db so first convert:
    SNR = 10 ** (SNRdb/10)
    ## Now use the definition: SNR = signal^2 / sigma^2
    sigma = np.sqrt(0.5 / SNR)
    #print(sigma)
    noise = prng.normal(0, sigma, length)
    sigmaActual = np.sqrt((np.sum(noise ** 2)) / length)
    noisyVector = vector + noise
    return noisyVector, sigma, sigmaActual