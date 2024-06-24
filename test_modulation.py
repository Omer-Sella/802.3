# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 13:51:25 2024

@author: Omer
"""

from modulationFunctions import modulatePAM4, pam4ClassifierToBits, pam4Slice, pam4SymbolsToBits, precoder, pam4ClassifierToBits
import numpy as np
from ieeeConstants import PAM4_LEVELS

def test_pam4Modulation():
    # From https://www.ieee802.org/3/bj/public/sep12/lusted_3bj_01_0912.pdf slide 25
    # Note that the presentation has some comments on underlined numbers, which is why I only used the first 40 bits == first 10 PAM4 symbols
    pmdLane0 = np.array([1,1,1,1,0,1,1,0,1,1,0,0,0,0,0,1,1,1,1,0, 0,1,0,1,1, 0,0,1,1,0, 0,1,0,1,1, 1,1,1,1,1, 1,1,1,0,1, 0,1,1,0,1, 1,0,11,1, 0,1,1,1,1, 0,0,1,1,1, 0,0,0,1,1, 0,0,1,0,0, 0,1,0,1,0, 0,1,0,0,0, 1,1,1,0,1])
    pmdLane0 = pmdLane0[0:40]
    pmdLan0Pam4SymbolsGrayCoded = np.array([2,2,1,3,2,0,0,1,2,3,1,1,3,1,3,1,1,2,2,2,2,3,3,1,3,2,1,2,1,2,3,1,2,0,1,3,1,0,1,1,0,3,0,2,3,3])
    
    pmdLane0Precoded = np.array([2,0,1,2,0,0,0,1,1,2,3,2,1,0,3,2,3,3,3,3,3,0,3,2,1,1,0,2,3,3,0,1,1,3,2,1,0,0,1,0,0,3,1,1,2,3])
    
    modulatedVector, pam4Symbols, pam4SymbolsPrecoded = modulatePAM4(pmdLane0, grayCoding = True, precoding = True, precoderInit = 2)
    assert (np.all(pmdLan0Pam4SymbolsGrayCoded[0:10] == pam4Symbols[0:10]))
    assert (np.all(pmdLane0Precoded[0:10] == pam4SymbolsPrecoded[0:10]))
    
    pmdLane1 = np.array([ 1,1,1,0,0, 0,0,0,0,0, 0,0,0,0,1, 0,1,1,0,1, 1,0,1,1,0, 0,0,0,1,0, 1,0,0,0,1, 0,0,0,0,1, 1,0,1,0,1, 0,1,0,1,1, 1,1,1,0,1, 1,0,1,1,0, 0,0,0,0,1, 1,1,1,0,0, 1,0,1,1,0, 0,1,1,0,0, 1,0,1,1,1, 1,1,1,1,1])
    pmdLane1 = pmdLane1[0:40]
    pmdLane1Pam4SymbolsGrayCoded = np.array([2,3,0,0,0,0,0,3,2,1,3,2,0,0,3,3,0,3,0,1,3,3,3,1,1,2,2,1,3,2,0,0,1,2,3,1,1,3,1,3,1,1,2,2,2,3])
    pmdLane1Precoded = np.array([2,1,3,1,3,1,3,0,2,3,0,2,2,2,1,2,2,1,3,2,1,2,3,2,3,3,3,2,1,1,3,1,0,2,1,0,1,2,3,0,1,0,2,0,2,3])
                                 
    modulatedVector, pam4Symbols, pam4SymbolsPrecoded = modulatePAM4(pmdLane1, grayCoding = True, precoding = True, precoderInit = 2)
    assert (np.all(pmdLane1Pam4SymbolsGrayCoded[0:10] == pam4Symbols[0:10]))
    assert (np.all(pmdLane1Precoded[0:10] == pam4SymbolsPrecoded[0:10]))
    
    pmdLane2 = np.array([  0,0,1,0,1, 1,1,0,1,0, 0,1,1,1,0, 1,0,1,0,1, 1,1,1,1,0, 0,0,0,0,0, 0,1,1,0,1, 0,0,0,0,0, 0,0,1,0,0, 0,0,1,1,0, 1,1,0,0,0, 1,1,0,0,1, 0,1,0,0,1, 1,1,1,1,1, 1,0,1,1,1, 1,1,0,1,1, 0,1,0,1,0, 1,1,0,0,0])
    pmdLane2 = pmdLane2[0:40]
    pmdLane2Pam4SymbolsGrayCoded = np.array([ 0,3,2,3,3,1,2,1,1,1,2,2,0,0,0,1,3,3,0,0,0,3,0,0,2,1,3,0,2,0,3,3,1,2,2,2,1,2,2,1,3,3,3,2,0,0])
    pmdLane2Precoded = np.array([0,3,3,0,3,2,0,1,0,1,1,1,3,1,3,2,1,2,2,2,2,1,0,0,2,3,0,0,2,2,1,2,3,3,3,3,2,0,2,3,0,3,0,2,2,0])
    
    modulatedVector, pam4Symbols, pam4SymbolsPrecoded = modulatePAM4(pmdLane2, grayCoding = True, precoding = True, precoderInit = 0)
    assert (np.all(pmdLane2Pam4SymbolsGrayCoded[0:10] == pam4Symbols[0:10]))
    assert (np.all(pmdLane2Precoded[0:10] == pam4SymbolsPrecoded[0:10]))
    
    pmdLane3 = np.array([0,0,1,0,1, 0,1,1,1,1, 0,0,1,0,1, 1,0,1,0,1, 1,1,1,1,0, 1,0,0,0,0, 0,1,1,0,1, 0,1,1,1,0, 1,1,1,1,1, 0,1,1,1,0, 0,0,1,1,0, 1,0,1,0,1, 0,0,0,1,0, 0,0,0,0,0, 0,0,1,1,1, 0,0,0,0,0, 0,0,0,1,0, 1,0,0,0,0])
    pmdLane3 = pmdLane3[0:40]
    pmdLane3Pam4SymbolsGrayCoded = np.array([ 0,3,3,2,2,0,3,2,1,1,2,2,1,0,0,1,3,3,2,3,2,2,3,1,2,0,1,3,3,3,3,0,3,0,0,0,1,2,0,0,0,0,3,3,0,0])
    pmdLane3Precoded = np.array([ 0,3,0,2,0,0,3,3,2,3,3,3,2,2,2,3,0,3,3,0,2,0,3,2,0,0,1,2,1,2,1,3,0,0,0,0,1,1,3,1,3,1,2,1,3,0])
    
    modulatedVector, pam4Symbols, pam4SymbolsPrecoded = modulatePAM4(pmdLane3, grayCoding = True, precoding = True, precoderInit = 0)
    assert (np.all(pmdLane3Pam4SymbolsGrayCoded[0:10] == pam4Symbols[0:10]))
    assert (np.all(pmdLane3Precoded[0:10] == pam4SymbolsPrecoded[0:10]))
    
    
def test_prbs9():
    from ieeeConstants import PRBS9Q_seed111111111_clause_120_5_11_2_a
    pass
    

def test_roundrtipGrayCoded():
    bitsTx = np.random.randint(0,2,1000)
    modulatedVector, pam4SymbolsTx, pam4SymbolsPrecoded = modulatePAM4(bitsTx, grayCoding = True, precoding = False)
    pam4SymbolsRx, errorAbsoluteValue, classifier = pam4Slice(modulatedVector)
    bitsRx = pam4SymbolsToBits(pam4SymbolsRx, grayCoded = True)
    assert (np.all(bitsTx == bitsRx))
    
def test_roundrtipNoGrayCoding():
    bitsTx = np.random.randint(0,2,1000)
    modulatedVector, pam4SymbolsTx, pam4SymbolsPrecoded = modulatePAM4(bitsTx, grayCoding = False, precoding = False)
    pam4SymbolsRx, errorAbsoluteValue, classifier = pam4Slice(modulatedVector)
    bitsRx = pam4SymbolsToBits(pam4SymbolsRx, grayCoded = False)
    assert (np.all(bitsTx == bitsRx))
    
def test_pam4ClassifierOverData():
    # Classifier is a 4 X len(modulatedVector) array such that: 
    # argmax over axis==0 should give the pam4 symbol
    # the sum over axis=0 (so vertically, i.e. the sum of 4 possibilities) == 1
    # all elements are >= 0
    bitsTx = np.random.randint(0,2,1000)
    modulatedVector, pam4SymbolsTx, pam4SymbolsPrecoded = modulatePAM4(bitsTx, grayCoding = False, precoding = False)
    pam4SymbolsRx, errorAbsoluteValue, classifier = pam4Slice(modulatedVector)
    assert(np.all(pam4SymbolsRx == np.argmax(classifier, axis = 0)))
    tolerance = 0.1
    assert(np.all( (np.sum(classifier, axis = 0) - 1) < tolerance) )
    assert (np.all(classifier >= 0))
    
def test_pam4ClassifierOverRandomNumbers(sampleSize = 1000):
    # Classifier is a 4 X len(modulatedVector) array such that: 
    # argmax over axis==0 should give the pam4 symbol
    # the sum over axis=0 (so vertically, i.e. the sum of 4 possibilities) == 1
    # all elements are >= 0
    
    modulatedVector = np.random.uniform(low = -3, high = 3, size = sampleSize)
    pam4SymbolsRx, errorAbsoluteValue, classifier = pam4Slice(modulatedVector)
    assert(np.all(pam4SymbolsRx == np.argmax(classifier, axis = 0)))
    tolerance = 0.1
    assert(np.all((np.sum(classifier, axis = 0) - 1) < tolerance))
    assert (np.all(classifier >= 0))

def test_pam4ClassifierToBits(sampleSize = 1):
    levels = np.array([])
    modulatedVector = np.random.choice(PAM4_LEVELS, size = sampleSize)
    pam4SymbolsRx, errorAbsoluteValue, classifier = pam4Slice(modulatedVector)
    # Safety - the purpose of the test is to create perfect conditions except for one index
    for i in range(sampleSize):
        assert(np.sum(classifier[:,i]) == 1)
        print(np.count_nonzero(classifier[:,i]))
        assert(np.count_nonzero(classifier[:,i]) == 1)
    
    for i in range(sampleSize):
        testIndex = np.random.randint(0, sampleSize)
        # This should give 0.5 for 0 and 0.5 for 1 
        classifier[testIndex,:] = np.array([0.25,0.25,0.25,0.25])
        bits, probabilities = pam4ClassifierToBits(classifier)
        assert(probabilities[2 * testIndex] == 0.5)
        assert(probabilities[2 * testIndex + 1] == 0.5)
        
    for i in range(sampleSize):
        testIndex = np.random.randint(0, sampleSize)
        # This should give 0.5 for 0 and 0.5 for 1 
        classifier[testIndex,:] = np.array([0.7,0.1,0.1,0.1])
        bits, probabilities = pam4ClassifierToBits(classifier)
        assert(probabilities[2 * testIndex] == 0.8)
        assert(probabilities[2 * testIndex + 1] == 0.8)
    


def test_precoder():
    """
    Precoder test according to the test sequences in parthasarathy_01_0911.pdf
    """
    testSequence = np.array([2, 2, 2, 2, 0, 3, 2, 0, 1, 3, 3, 0, 0, 0, 0, 2, 3, 0, 3]) # Taken from the presentation as is
    desiredOutput = np.array([0, 2, 0, 2, 2, 1, 1, 3, 2, 1, 2, 2, 2, 2, 2, 0, 3, 1, 2]) # Taken from the presentation as is
    precoderOutput = precoder(testSequence)
    assert(np.all(precoderOutput == desiredOutput))
  
