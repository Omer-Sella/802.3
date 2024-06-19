# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 13:05:28 2023

@author: Megatron
"""
import os
import numpy as np
#from numba import jit, int32, float32, types, typed, boolean, float64, int64
#import math

projectDir = os.environ.get('8023')
if projectDir == None:
    import pathlib
    projectDir = pathlib.Path(__file__).parent.absolute()

import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, projectDir)
#import io


def additiveWhiteGaussianNoise(vector, length, SNRdb, prng, d = 0.5):
    
    """
    The problem is this: I am passing SNRdb as a number, and wish to extract sigma.
    However, the formulas for SNR for PAM2 and for SNR for PAM4 are different.
    In PAM2 SNR == EbN0 (Energy pre bit / N0), and in PAM4 the average energy per bit (!) is different (continue in Proakis page 268)
    """
    ## The input SNR is in db so first convert:
    SNR = 10 ** (SNRdb/10)
    ## Now use the definition: SNR = signal^2 / sigma^2
    # Consult https://www.ieee802.org/3/ap/public/jul04/liu_01_0704.pdf 
    # and https://grouper.ieee.org/groups/802/3/bj/public/nov11/moore_01a_1111.pdf
    # and https://uk.mathworks.com/help/comm/ug/bit-error-rate-analysis-techniques.html#a1044894202
    # and https://www.dsprelated.com/showarticle/168.php
    # and https://uk.mathworks.com/help/comm/ref/berawgn.html
    # and https://uk.mathworks.com/help/comm/ug/analytical-expressions-used-in-berawgn-function-and-bit-error-rate-analysis-app.html
    
    
    sigma = np.sqrt(d / SNR)
    noise = prng.normal(0, sigma, length)
    sigmaActual = np.sqrt((np.sum(noise ** 2)) / length)
    noisyVector = vector + noise
    return noisyVector, sigma, sigmaActual

