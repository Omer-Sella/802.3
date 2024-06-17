# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 13:51:25 2024

@author: Omer
"""

from modulationFunctions import modulatePAM4
import numpy as np

def test_pam4Modulation():
    # From https://www.ieee802.org/3/bj/public/sep12/lusted_3bj_01_0912.pdf slide 25
    pmdLane0 = np.array([1,1,1,1,0,1,1,0,1,1,0,0,0,0,0,1,1,1,1,0, 0,1,0,1,1, 0,0,1,1,0, 0,1,0,1,1, 1,1,1,1,1, 1,1,1,0,1, 0,1,1,0,1, 1,0,1,1,1, 0,1,1,1,1, 0,0,1,1,1, 0,0,0,1,1, 0,0,1,0,0, 0,1,0,1,0, 0,1,0,0,0, 1,1,1,0,1])
    pmdLan0Pam4SymbolsGrayCoded = np.array([2,2,1,3,2,0,0,1,2,3,1,1,3,1,3,1,1,2,2,2,2,3,3,1,3,2,1,2,1,2,3,1,2,0,1,3,1,0,1,1,0,3,0,2,3,3])
    
    pmdLane0Precoded = np.array([2,0,1,2,0,0,0,1,1,2,3,2,1,0,3,2,3,3,3,3,3,0,3,2,1,1,0,2,3,3,0,1,1,3,2,1,0,0,1,0,0,3,1,1,2,3])
    
    modulatedVector, pam4Symbols, pam4SymbolsPrecoded = modulatePAM4(pmdLane0, grayCoding = True, precoding = True)
    #assert (np.all(pmdLan0Pam4SymbolsGrayCoded == pam4Symbols))
    #assert (np.all(pmdLane0Precoded == pam4SymbolsPrecoded))
    
    pmdLane1 = np.array([ 1,1,1,0,0, 0,0,0,0,0, 0,0,0,0,1, 0,1,1,0,1, 1,0,1,1,0, 0,0,0,1,0, 1,0,0,0,1, 0,0,0,0,1, 1,0,1,0,1, 0,1,0,1,1, 1,1,1,0,1, 1,0,1,1,0, 0,0,0,0,1, 1,1,1,0,0, 1,0,1,1,0, 0,1,1,0,0, 1,0,1,1,1, 1,1,1,1,1])
    pmdLane1Pam4SymbolsGrayCoded = np.array([2,3,0,0,0,0,0,3,2,1,3,2,0,0,3,3,0,3,0,1,3,3,3,1,1,2,2,1,3,2,0,0,1,2,3,1,1,3,1,3,1,1,2,2,2,3])
    pmdLane1Precoded = np.array([2,3,0,0,0,0,0,3,2,1,3,2,0,0,3,3,0,3,0,1,3,3,3,1,1,2,2,1,3,2,0,0,1,2,3,1,1,3,1,3,1,1,2,2,2,3])
    
    modulatedVector, pam4Symbols, pam4SymbolsPrecoded = modulatePAM4(pmdLane1, grayCoding = True, precoding = True)
    assert (np.all(pmdLane1Pam4SymbolsGrayCoded == pam4Symbols))
    assert (np.all(pmdLane1Precoded == pam4SymbolsPrecoded))
    
    pmdLane2 = np.array([  0,0,1,0,1, 1,1,0,1,0, 0,1,1,1,0, 1,0,1,0,1, 1,1,1,1,0, 0,0,0,0,0, 0,1,1,0,1, 0,0,0,0,0, 0,0,1,0,0, 0,0,1,1,0, 1,1,0,0,0, 1,1,0,0,1, 0,1,0,0,1, 1,1,1,1,1, 1,0,1,1,1, 1,1,0,1,1, 0,1,0,1,0, 1,1,0,0,0])
    pmdLane2Pam4SymbolsGrayCoded = np.array([ 0,3,2,3,3,1,2,1,1,1,2,2,0,0,0,1,3,3,0,0,0,3,0,0,2,1,3,0,2,0,3,3,1,2,2,2,1,2,2,1,3,3,3,2,0,0])
    pmdLane2Precoded = np.array([0,3,3,0,3,2,0,1,0,1,1,1,3,1,3,2,1,2,2,2,2,1,0,0,2,3,0,0,2,2,1,2,3,3,3,3,2,0,2,3,0,3,0,2,2,0])
    
    modulatedVector, pam4Symbols, pam4SymbolsPrecoded = modulatePAM4(pmdLane2, grayCoding = True, precoding = True)
    assert (np.all(pmdLane2Pam4SymbolsGrayCoded == pam4Symbols))
    assert (np.all(pmdLane2Precoded == pam4SymbolsPrecoded))
    
    pmdLane3 = np.array([0,0,1,0,1, 0,1,1,1,1, 0,0,1,0,1, 1,0,1,0,1, 1,1,1,1,0, 1,0,0,0,0, 0,1,1,0,1, 0,1,1,1,0, 1,1,1,1,1, 0,1,1,1,0, 0,0,1,1,0, 1,0,1,0,1, 0,0,0,1,0, 0,0,0,0,0, 0,0,1,1,1, 0,0,0,0,0, 0,0,0,1,0, 1,0,0,0,0])
    pmdLane3Pam4SymbolsGrayCoded = np.array([ 0,3,3,2,2,0,3,2,1,1,2,2,1,0,0,1,3,3,2,3,2,2,3,1,2,0,1,3,3,3,3,0,3,0,0,0,1,2,0,0,0,0,3,3,0,0])
    pmdLane3Precoded = np.array([ 0,3,0,2,0,0,3,3,2,3,3,3,2,2,2,3,0,3,3,0,2,0,3,2,0,0,1,2,1,2,1,3,0,0,0,0,1,1,3,1,3,1,2,1,3,0])
    
    modulatedVector, pam4Symbols, pam4SymbolsPrecoded = modulatePAM4(pmdLane3, grayCoding = True, precoding = True)
    assert (np.all(pmdLane3Pam4SymbolsGrayCoded == pam4Symbols))
    assert (np.all(pmdLane3Precoded == pam4SymbolsPrecoded))
    
    
def test_prbs9():
    from ieeeConstrants import PRBS9Q_seed111111111_clause_120_5_11_2_a
    pass
    
    
    