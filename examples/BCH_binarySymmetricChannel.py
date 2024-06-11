# -*- coding: utf-8 -*-
"""
Binary Symmetric channel decoding

Created on Tue Jun 11 10:20:30 2024

@author: Omer
"""

import os, sys
reedSolomonProjectDir = os.environ.get('REEDSOLOMON')
if reedSolomonProjectDir == None: 
     reedSolomonProjectDir = "c:/users/omer/reedSolomon/reedSolomon/"
sys.path.insert(0, reedSolomonProjectDir)
import os, sys
ieeeProjectDir = os.environ.get('IEEE8023')
if ieeeProjectDir == None: 
     ieeeProjectDir = "c:/users/omer/802.3/"
sys.path.insert(0, ieeeProjectDir)
import numpy as np
from arithmetic import generateExponentAndLogTables, polynomial, gf128, binaryFieldElement as gf2
from bchDecoder import bchDecoder
from ieee8023dj_d0p1 import bchEncoder
from graphics import plotSNRvsBER
import modulationFunctions as mf
import channelFunctions as cf
from ieeeConstants import snrBaseline, berPam2
localPrng = np.random.RandomState(seed = 8023)

eD, _ =  generateExponentAndLogTables()
generatorMatrix = np.load("c:/users/omer/802.3/bchMatrixEncoder.npy")


bscBerStats = []
numberOfRepeats = 10
for p in berPam2:  
    error = 0
    for i in range(numberOfRepeats):
        binaryData = np.random.randint(0,2,110)
        encodedBinaryData = bchEncoder(binaryData, G = generatorMatrix)
        errorVector = np.random.binomial(1,p,126)
        encodedBinaryData = (encodedBinaryData + error) %2
        correctedVector, correctionVector, errorLocatorX = bchDecoder( receivedBinaryVecotor = encodedBinaryData,
                                                                      exponentDictionary = eD,
                                                                      numberOfPowers = 16,
                                                                      codewordLengthMaximal = 127)
        
        error = error + np.sum(correctedVector != encodedBinaryData)
        print(np.sum(correctedVector != encodedBinaryData))
    bscBerStats.append(error)
print(bscBerStats)
print(bscBerStats / (numberOfRepeats * 126))