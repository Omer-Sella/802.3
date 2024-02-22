# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 10:39:30 2024

@author: Omer
"""
import os
import sys
projectDir = os.environ.get('IEEE8032DJ')
if projectDir == None:
     projectDir = "D:/802.3/"
sys.path.insert(1, projectDir)
import numpy as np
from ieee8023dj_d0p1 import *
from modulationFunctions import *
from channelFunctions import *
from graphics import *


def simpleAWGNComparison(dataLength = 500, seed = 7134066, snrAxis = [2,3,4,5,6]):
    localRandom = np.random.RandomState(seed)    
    #Generate random data
    randomData = localRandom.randint(0, 1, size = 2 * dataLength)
    
    #Modulate
    pam4ModulatedNoGreyCoding = modulatePAM4(randomData, greyCoded = False)
    pam4ModulatedGreyCoded = modulatePAM4(randomData, greyCoded = True)
    nrzModulated = modulatePAM2(randomData)
    
    berPAM4Coded = np.ones(len(snrAxis))
    berPAM4NoCoding = np.ones(len(snrAxis))
    berNrz = np.ones(len(snrAxis))
    for i in range(len(snrAxis)):
        
        #Add noise
        SNRdb = snrAxis[i]
        pam4ModulatedNoGreyCodingNoisy, _ , _ = additiveWhiteGaussianNoise(pam4ModulatedNoGreyCoding, pam4ModulatedNoGreyCoding.shape[0], SNRdb, localRandom)
        pam4ModulatedGreyCodedNoisy, _ , _ = additiveWhiteGaussianNoise(pam4ModulatedGreyCoded, pam4ModulatedNoGreyCoding.shape[0], SNRdb, localRandom)
        nrzModulatedNoisy, _ , _ = additiveWhiteGaussianNoise(nrzModulated, nrzModulated.shape[0], SNRdb, localRandom)
        
        #Slice
        pam4ModulatedNoGreyCodingNoisySliced, _ = pam4Slicer(pam4ModulatedNoGreyCodingNoisy, greyCoded = False)
        pam4ModulatedGreyCodedNoisySliced, _ = pam4Slicer(pam4ModulatedGreyCodedNoisy, greyCoded = True)
        nrzModulatedNoisySliced, _ = slicer(nrzModulatedNoisy)
        
        #Calculate BER
        berPAM4NoCoding[i] = np.count(pam4ModulatedNoGreyCodingNoisySliced != randomData) / (2 * dataLength)
        berPAM4Coded[i] = np.count(pam4ModulatedGreyCodedNoisySliced != randomData) / (2 * dataLength)
        berNrz[i] = np.count(nrzModulatedNoisySliced != randomData) / (2 * dataLength)
    
    return berNrz, berPAM4Coded, berPAM4NoCoding

