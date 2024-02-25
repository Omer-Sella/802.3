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


def simpleAWGNComparison(dataLength = 500, seed = 7134066, snrAxis = [2,3,4,5,6,7,8,9,10]):
    localRandom = np.random.RandomState(seed)    
    #Generate random data
    randomData = localRandom.randint(0, 1, size = 2 * dataLength)
    #print("Random data\n")
    #print(randomData)
    #Modulate
    pam4ModulatedNoGreyCoding = modulatePAM4(randomData, greyCoded = False)
    #print("No grey coding")
    #print(pam4ModulatedNoGreyCoding)
    pam4ModulatedGreyCoded = modulatePAM4(randomData, greyCoded = True)
    #print("With grey coding")
    #print(pam4ModulatedGreyCoded)
    nrzModulated = modulatePAM2(randomData)
    
    berPAM4Coded = np.ones(len(snrAxis))
    berPAM4NoCoding = np.ones(len(snrAxis))
    berNrz = np.ones(len(snrAxis))
    for i in range(len(snrAxis)):
        
        #Add noise
        SNRdb = snrAxis[i]
        zro = np.zeros(pam4ModulatedGreyCoded.shape[0])
        noise,_ ,_  = additiveWhiteGaussianNoise(zro, zro.shape[0], SNRdb, localRandom)
        pam4ModulatedNoGreyCodingNoisy = pam4ModulatedNoGreyCoding + noise
        
        #print(pam4ModulatedNoGreyCodingNoisy)
        pam4ModulatedGreyCodedNoisy = pam4ModulatedGreyCoded + noise
        #assert(np.all(pam4ModulatedNoGreyCodingNoisy == pam4ModulatedGreyCodedNoisy))
        #assert(np.all(pam4ModulatedNoGreyCodingNoisy == pam4ModulatedGreyCodedNoisy))
        #pam4ModulatedNoGreyCodingNoisy, _ , _ = additiveWhiteGaussianNoise(pam4ModulatedNoGreyCoding, pam4ModulatedNoGreyCoding.shape[0], SNRdb, localRandom)
        #pam4ModulatedGreyCodedNoisy, _ , _ = additiveWhiteGaussianNoise(pam4ModulatedGreyCoded, pam4ModulatedGreyCoded.shape[0], SNRdb, localRandom)
        nrzModulatedNoisy, _ , _ = additiveWhiteGaussianNoise(nrzModulated, nrzModulated.shape[0], SNRdb, localRandom)
        
        #Slice
        pam4ModulatedNoGreyCodingNoisySliced, debug = pam4Slicer(pam4ModulatedNoGreyCodingNoisy, greyCoded = False)
        #print(debug)
        #print(pam4ModulatedNoGreyCodingNoisySliced)
        pam4ModulatedGreyCodedNoisySliced, debug = pam4Slicer(pam4ModulatedGreyCodedNoisy, greyCoded = True)
        #print(debug)
        #print(pam4ModulatedGreyCodedNoisySliced)
        nrzModulatedNoisySliced = slicer(nrzModulatedNoisy)
        
        
        #Calculate BER
        berPAM4NoCoding[i] = np.count_nonzero(pam4ModulatedNoGreyCodingNoisySliced != randomData) / (2 * dataLength)
        berPAM4Coded[i] = np.count_nonzero(pam4ModulatedGreyCodedNoisySliced != randomData) / (2 * dataLength)
        berNrz[i] = np.count_nonzero(nrzModulatedNoisySliced != randomData) / (2 * dataLength)
    
    return berNrz, berPAM4Coded, berPAM4NoCoding

