# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 14:26:43 2024

@author: Omer
"""

import numpy as np
from ieeeConstants import *
from ieee8023dj_d0p1 import bchEncoder as ieeeBchEncoder
import os
import sys
projectDir = os.environ.get('IEEE8032DJ')
#You don't have to define an environment variable, but then you have to give a path to the project here
if projectDir == None: 
     projectDir = "c:/users/omer/802.3/"
sys.path.insert(1, projectDir)
import numpy as np
from ieee8023dj_d0p1 import bchEncoder as ieeeBchEncoder, encodeUsingMatrixOnly

def test_bchEncoder(numberOfAttempts = 1000):
    # This function tests / verifies / validates that the encoder function implementation is identical to a matrix generated indepedently by a different bch encoder in a project named reedSolomon
    pathToGeneratorMatrix = projectDir + "/bchMatrixEncoder.npy"
    generatorMatrix = np.load(pathToGeneratorMatrix)
    for i in range(numberOfAttempts):
        data = np.random.randint(0,2,110)
        encodedData1 = ieeeBchEncoder(data, G = None)
        encodedData2 = ieeeBchEncoder(data, G = generatorMatrix)
        assert(np.all(encodedData1 == encodedData2))

def test_bchEncoder_zero_data():
    zeroData = np.zeros(110, dtype = IEEE_8023_INT_DATA_TYPE)
    #tv0[109] = 1
    zeroDataEncoded = ieeeBchEncoder(zeroData)
    assert (np.all(zeroDataEncoded == 0))
            