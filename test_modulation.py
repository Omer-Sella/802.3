# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 13:51:25 2024

@author: Omer
"""

from modulationFunctions import modulatePAM4, pam4Slice, pam4SymbolsToBits
import numpy as np

def test_pam4Modulation():
    # From https://www.ieee802.org/3/bj/public/sep12/lusted_3bj_01_0912.pdf slide 25
    # Note that the presentation has some comments on underlined numbers, which is why I only used the first 40 bits == first 10 PAM4 symbols
    pmdLane0 = np.array([1,1,1,1,0,1,1,0,1,1,0,0,0,0,0,1,1,1,1,0, 0,1,0,1,1, 0,0,1,1,0, 0,1,0,1,1, 1,1,1,1,1, 1,1,1,0,1, 0,1,1,0,1, 1,0,11,1, 0,1,1,1,1, 0,0,1,1,1, 0,0,0,1,1, 0,0,1,0,0, 0,1,0,1,0, 0,1,0,0,0, 1,1,1,0,1])
    pmdLane0 = pmdLane0[0:40]
    pmdLan0Pam4SymbolsGrayCoded = np.array([2,2,1,3,2,0,0,1,2,3,1,1,3,1,3,1,1,2,2,2,2,3,3,1,3,2,1,2,1,2,3,1,2,0,1,3,1,0,1,1,0,3,0,2,3,3])
    
    pmdLane0Precoded = np.array([2,0,1,2,0,0,0,1,1,2,3,2,1,0,3,2,3,3,3,3,3,0,3,2,1,1,0,2,3,3,0,1,1,3,2,1,0,0,1,0,0,3,1,1,2,3])
    
    modulatedVector, pam4Symbols, pam4SymbolsPrecoded = modulatePAM4(pmdLane0, grayCoding = True, precoding = True)
    assert (np.all(pmdLan0Pam4SymbolsGrayCoded[0:10] == pam4Symbols[0:10]))
    #assert (np.all(pmdLane0Precoded[0:10] == pam4SymbolsPrecoded[0:10]))
    
    pmdLane1 = np.array([ 1,1,1,0,0, 0,0,0,0,0, 0,0,0,0,1, 0,1,1,0,1, 1,0,1,1,0, 0,0,0,1,0, 1,0,0,0,1, 0,0,0,0,1, 1,0,1,0,1, 0,1,0,1,1, 1,1,1,0,1, 1,0,1,1,0, 0,0,0,0,1, 1,1,1,0,0, 1,0,1,1,0, 0,1,1,0,0, 1,0,1,1,1, 1,1,1,1,1])
    pmdLane1 = pmdLane1[0:40]
    pmdLane1Pam4SymbolsGrayCoded = np.array([2,3,0,0,0,0,0,3,2,1,3,2,0,0,3,3,0,3,0,1,3,3,3,1,1,2,2,1,3,2,0,0,1,2,3,1,1,3,1,3,1,1,2,2,2,3])
    pmdLane1Precoded = np.array([2,3,0,0,0,0,0,3,2,1,3,2,0,0,3,3,0,3,0,1,3,3,3,1,1,2,2,1,3,2,0,0,1,2,3,1,1,3,1,3,1,1,2,2,2,3])
    
    modulatedVector, pam4Symbols, pam4SymbolsPrecoded = modulatePAM4(pmdLane1, grayCoding = True, precoding = True)
    assert (np.all(pmdLane1Pam4SymbolsGrayCoded[0:10] == pam4Symbols[0:10]))
    #assert (np.all(pmdLane1Precoded[0:10] == pam4SymbolsPrecoded[0:10]))
    
    pmdLane2 = np.array([  0,0,1,0,1, 1,1,0,1,0, 0,1,1,1,0, 1,0,1,0,1, 1,1,1,1,0, 0,0,0,0,0, 0,1,1,0,1, 0,0,0,0,0, 0,0,1,0,0, 0,0,1,1,0, 1,1,0,0,0, 1,1,0,0,1, 0,1,0,0,1, 1,1,1,1,1, 1,0,1,1,1, 1,1,0,1,1, 0,1,0,1,0, 1,1,0,0,0])
    pmdLane2 = pmdLane2[0:40]
    pmdLane2Pam4SymbolsGrayCoded = np.array([ 0,3,2,3,3,1,2,1,1,1,2,2,0,0,0,1,3,3,0,0,0,3,0,0,2,1,3,0,2,0,3,3,1,2,2,2,1,2,2,1,3,3,3,2,0,0])
    pmdLane2Precoded = np.array([0,3,3,0,3,2,0,1,0,1,1,1,3,1,3,2,1,2,2,2,2,1,0,0,2,3,0,0,2,2,1,2,3,3,3,3,2,0,2,3,0,3,0,2,2,0])
    
    modulatedVector, pam4Symbols, pam4SymbolsPrecoded = modulatePAM4(pmdLane2, grayCoding = True, precoding = True)
    assert (np.all(pmdLane2Pam4SymbolsGrayCoded[0:10] == pam4Symbols[0:10]))
    #assert (np.all(pmdLane2Precoded[0:10] == pam4SymbolsPrecoded[0:10]))
    
    pmdLane3 = np.array([0,0,1,0,1, 0,1,1,1,1, 0,0,1,0,1, 1,0,1,0,1, 1,1,1,1,0, 1,0,0,0,0, 0,1,1,0,1, 0,1,1,1,0, 1,1,1,1,1, 0,1,1,1,0, 0,0,1,1,0, 1,0,1,0,1, 0,0,0,1,0, 0,0,0,0,0, 0,0,1,1,1, 0,0,0,0,0, 0,0,0,1,0, 1,0,0,0,0])
    pmdLane3 = pmdLane3[0:40]
    pmdLane3Pam4SymbolsGrayCoded = np.array([ 0,3,3,2,2,0,3,2,1,1,2,2,1,0,0,1,3,3,2,3,2,2,3,1,2,0,1,3,3,3,3,0,3,0,0,0,1,2,0,0,0,0,3,3,0,0])
    pmdLane3Precoded = np.array([ 0,3,0,2,0,0,3,3,2,3,3,3,2,2,2,3,0,3,3,0,2,0,3,2,0,0,1,2,1,2,1,3,0,0,0,0,1,1,3,1,3,1,2,1,3,0])
    
    modulatedVector, pam4Symbols, pam4SymbolsPrecoded = modulatePAM4(pmdLane3, grayCoding = True, precoding = True)
    assert (np.all(pmdLane3Pam4SymbolsGrayCoded[0:10] == pam4Symbols[0:10]))
    #assert (np.all(pmdLane3Precoded[0:10] == pam4SymbolsPrecoded[0:10]))
    
    
def test_prbs9():
    from ieeeConstants import PRBS9Q_seed111111111_clause_120_5_11_2_a
    pass
    

def test_roundrtipGrayCoded():
    bitsTx = np.random.randint(0,2,1000)
    modulatedVector, pam4SymbolsTx, pam4SymbolsPrecoded = modulatePAM4(bitsTx, grayCoding = True, precoding = False)
    pam4SymbolsRx, errorAbsoluteValue = pam4Slice(modulatedVector)
    bitsRx = pam4SymbolsToBits(pam4SymbolsRx, grayCoded = True)
    assert (np.all(bitsTx == bitsRx))
    
def test_roundrtipNoGrayCoding():
    bitsTx = np.random.randint(0,2,1000)
    modulatedVector, pam4SymbolsTx, pam4SymbolsPrecoded = modulatePAM4(bitsTx, grayCoding = False, precoding = False)
    pam4SymbolsRx, errorAbsoluteValue = pam4Slice(modulatedVector)
    bitsRx = pam4SymbolsToBits(pam4SymbolsRx, grayCoded = False)
    assert (np.all(bitsTx == bitsRx))
    
    