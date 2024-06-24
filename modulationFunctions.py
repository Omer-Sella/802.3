# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 13:03:02 2023

@author: Megatron
"""
import os
import sys
#from numba import jit, int32, float32, types, typed, boolean, float64, int64
#import math
projectDir = os.environ.get('IEEE8032DJ')
if projectDir == None:
    import pathlib
    projectDir = pathlib.Path(__file__).parent.absolute()
sys.path.insert(1, projectDir)
import numpy as np
import time
import ieeeConstants
from ieeeConstants import IEEE_8023_INT_DATA_TYPE, IEEE_8023_DECIMAL_DATA_TYPE, PAM4_GRAYCODED, PAM4_LEVEL_HIGH, PAM4_LEVEL_LOW, PAM4_LEVEL_MID_HIGH, PAM4_LEVEL_MID_LOW, PAM4_NOGRAYCODING, PAM4_GRAYCODED



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

def precoder(vector, init = 'default'):
    # 1/(1+D) modulu 4 precoder - need referece
    vectorPrecoded = np.zeros(len(vector), dtype = IEEE_8023_INT_DATA_TYPE)
    if init == 'default':
        vectorPrecoded[0] = (vector[0] - vector[1])%4
    else:
        assert(init == 0 or init == 1 or init == 2 or init == 3)
        vectorPrecoded[0] = init
    for i in range(1,len(vector)):
        vectorPrecoded[i] = (vector[i] - vectorPrecoded[i-1]) %4
    return vectorPrecoded
    
    

def modulatePAM4(vector, grayCoding = True, precoding = False, precoderInit = 'default', 
                 levels = [ieeeConstants.PAM4_LEVEL_LOW,
                           ieeeConstants.PAM4_LEVEL_MID_LOW, 
                           ieeeConstants.PAM4_LEVEL_MID_HIGH, 
                           ieeeConstants.PAM4_LEVEL_HIGH]):
    # 177.4.7.1 of draft 1.0 refers to 120.5.7.1 
    # 120.5.7.1 Gray mapping for PAM4 encoded lanes
    # For output lanes encoded as PAM4 (for 200GBASE-R, where the number of output lanes is 4, or for 
    # 400GBASE-R, where the number of output lanes is 4 or 8), the PMA transmit process shall map consecutive 
    # pairs of bits {A, B}, where A is the bit arriving first, to a Gray-coded symbol as follows:
    # {0, 0} maps to 0,
    # {0, 1} maps to 1,
    # {1, 1} maps to 2, and
    # {1, 0} maps to 3
    #vector is assumed to be a binary vector of even length
    # note that the following line means the stream goes [MSB0,LSB0,MSB1,LSB1 ...]
    if grayCoding:
        pam4Symbols = 2 * np.array(vector[0::2])  +  ((np.array(vector[1::2]) + np.array(vector[0::2])) %2)
    else:
        pam4Symbols = 2 * vector[0::2]  +  vector[1::2]
        
    if precoding:
        pam4SymbolsPrecoded = precoder(pam4Symbols, init = precoderInit)
    else:
        pam4SymbolsPrecoded = pam4Symbols
    
    modulatedVector = np.zeros(vector.shape[0] // 2, dtype = IEEE_8023_DECIMAL_DATA_TYPE)
    modulatedVector[np.where(pam4SymbolsPrecoded == 0)] = levels[0] #PAM4_LEVEL_LOW
    modulatedVector[np.where(pam4SymbolsPrecoded == 1)] = levels[1] #PAM4_LEVEL_MID_LOW
    modulatedVector[np.where(pam4SymbolsPrecoded == 2)] = levels[2] #PAM4_LEVEL_MID_HIGH
    modulatedVector[np.where(pam4SymbolsPrecoded == 3)] = levels[3] #PAM4_LEVEL_HIGH
    
    return modulatedVector, pam4Symbols, pam4SymbolsPrecoded

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

def pam4Quantize(signal, effectiveNumberOfBits = 6, high = ieeeConstants.PAM4_LEVEL_HIGH, low = ieeeConstants.PAM4_LEVEL_LOW):
    q = 1/(2 ** effectiveNumberOfBits)
    quantizedSignal = q * np.floor(signal/q)
    return quantizedSignal
    
def pam4Slice(signal, reference = np.array([PAM4_LEVEL_LOW, PAM4_LEVEL_MID_LOW, PAM4_LEVEL_MID_HIGH, PAM4_LEVEL_HIGH])):
    """
    Input signal of type np array of real numbers assumed to be a 1 dimensional column vector - no safety
    """
    referenceRepeated = np.tile(reference[:,np.newaxis], (1,len(signal)))
    signalRepeated = np.tile(signal, (len(reference),1))
    g = lambda x: np.exp(-1*np.abs(x))
    classifier =  g(signalRepeated - referenceRepeated) / np.sum( g(signalRepeated - referenceRepeated), axis = 0)
    # the result of np.argmin(np.abs(signalRepeated - referenceRepeated), axis = 0) is a vector with values 0,1,2,3 depending on the ordering of the reference levels - no safety, the user needs to think how they order the reference so that the PAM4 slicing will be correct
    pam4Symbols = np.argmin(np.abs(signalRepeated - referenceRepeated), axis = 0)
    errorAbsoluteValue = np.min(np.abs(signalRepeated - referenceRepeated), axis = 0)
    return pam4Symbols, errorAbsoluteValue, classifier

def pam4SymbolsToBits(pam4Symbols, grayCoded = True):
    bits = np.zeros(len(pam4Symbols) * 2, dtype = IEEE_8023_INT_DATA_TYPE)
    if grayCoded:
        for i in range(len(pam4Symbols)):
            bits[2 * i] = 1 * (pam4Symbols[i] > 1)
            bits[2 * i + 1] = 1 * (pam4Symbols[i] == 1 or pam4Symbols[i] == 2)
    else:
        for i in range(len(pam4Symbols)):
            bits[2 * i] = pam4Symbols[i] // 2
            bits[2 * i + 1] = pam4Symbols[i] % 2
    return bits
        

def pam4ClassifierToBits(classifier, grayCoded = True):
    """
    The idea is simple - we multiply the classifier value (probability) by the corresponding bit pair and sum.
    There is no need to use pam4SymbolsToBits, because the classifier represents probabilities for ALL pam4 symbols.
    In fact, since the is no control logic (if,else) or list comprehension here, this should be faster than pam4SymbolsToBits.
    """
    
    if grayCoded:
        symbolsToBits = PAM4_GRAYCODED
    else:
        symbolsToBits = PAM4_NOGRAYCODING
    
    msbRepeated = np.tile(symbolsToBits[:,0], (1,classifier.shape[1]))
    lsbRepeated = np.tile(symbolsToBits[:,1], (1,classifier.shape[1]))
    
    probabilityOfReceivingOneMsb = np.sum(msbRepeated * classifier, axis = 0)
    probabilityOfReceivingOneLsb = np.sum(lsbRepeated * classifier, axis = 0)
    # each vector is a row vector, so stack them MSB on top of LSB, and then flatten it using Fortran ordering, i.e.: column first, so stack[0,0], stack[0,1], stack[1,0], stack[1,1] ...
    probabilityOfReceivingOne = np.vstack((probabilityOfReceivingOneMsb, probabilityOfReceivingOneLsb)).ravel(order = 'F')
    
    bits = pam4SymbolsToBits(np.argmax(clasifier, axis = 0), grayCoded = grayCoded)
    probabilities = np.where(bits == 1, [probabilityOfReceivingOne, 1 - probabilityOfReceivingOne])
    
    return bits, probabilities
