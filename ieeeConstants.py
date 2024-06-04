# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 20:26:57 2024

@author: Omer
"""

import scipy.io
import sys
import os
projectDir = os.environ.get('IEEE8023DJ')
if projectDir == None:
     projectDir = "c:/users/omer/802.3/"
sys.path.insert(0, projectDir)

import numpy as np
IEEE_8023_DATA_TYPE = np.int64
IEEE_8023_INT_DATA_TYPE = np.int64
IEEE_8023_DECIMAL_DATA_TYPE = np.float64
IEEE_8023_DATA_TYPE = np.int64

# Omer Sella: seeds can be integers between 0 and 2**31 - 1
IEEE_8023_MAX_SEED = 2**31 - 1

#NUMBA_INT = int64
#NUMBA_FLOAT = float64
#NUMBA_BOOL = boolean

PAM4_LEVEL_LOW = -1
PAM4_LEVEL_MID_LOW = -(1/3)
PAM4_LEVEL_MID_HIGH = 1/3
PAM4_LEVEL_HIGH = 1

G = np.matrix([
    [1,0,0,1,0,1,0,0],
     [0,1,0,0,1,0,1,0],
 [0,0,1,0,0,1,0,1],
 [1,1,0,0,1,0,1,1],
 [1,0,1,1,1,1,0,0],
 [0,1,0,1,1,1,1,0],
 [0,0,1,0,1,1,1,1],
 [1,1,0,0,1,1,1,0],
 [0,1,1,0,0,1,1,1],
 [1,1,1,0,1,0,1,0],
 [0,1,1,1,0,1,0,1],
 [1,1,1,0,0,0,1,1],
 [1,0,1,0,1,0,0,0],
 [0,1,0,1,0,1,0,0],
 [0,0,1,0,1,0,1,0],
 [0,0,0,1,0,1,0,1],
 [1,1,0,1,0,0,1,1],
 [1,0,1,1,0,0,0,0],
 [0,1,0,1,1,0,0,0],
 [0,0,1,0,1,1,0,0],
 [0,0,0,1,0,1,1,0],
 [0,0,0,0,1,0,1,1],
 [1,1,0,1,1,1,0,0],
 [0,1,1,0,1,1,1,0],
 [0,0,1,1,0,1,1,1],
 [1,1,0,0,0,0,1,0],
 [0,1,1,0,0,0,0,1],
 [1,1,1,0,1,0,0,1],
 [1,0,1,0,1,1,0,1],
 [1,0,0,0,1,1,1,1],
 [1,0,0,1,1,1,1,0],
 [0,1,0,0,1,1,1,1],
 [1,1,1,1,1,1,1,0],
 [0,1,1,1,1,1,1,1],
 [1,1,1,0,0,1,1,0],
 [0,1,1,1,0,0,1,1],
 [1,1,1,0,0,0,0,0],
 [0,1,1,1,0,0,0,0],
 [0,0,1,1,1,0,0,0],
 [0,0,0,1,1,1,0,0],
 [0,0,0,0,1,1,1,0],
 [0,0,0,0,0,1,1,1],
[1,1,0,1,1,0,1,0],
[0,1,1,0,1,1,0,1],
[1,1,1,0,1,1,1,1],
[1,0,1,0,1,1,1,0],
[0,1,0,1,0,1,1,1],
[1,1,1,1,0,0,1,0],
[0,1,1,1,1,0,0,1],
[1,1,1,0,0,1,0,1],
[1,0,1,0,1,0,1,1],
[1,0,0,0,1,1,0,0],
[0,1,0,0,0,1,1,0],
[0,0,1,0,0,0,1,1],
[1,1,0,0,1,0,0,0],
[0,1,1,0,0,1,0,0],
[0,0,1,1,0,0,1,0],
[0,0,0,1,1,0,0,1],
[1,1,0,1,0,1,0,1],
[1,0,1,1,0,0,1,1]])

tv1_tp4 = np.array([0,1,0,1,1,0,0,0,1,1,0,0,1,1,0,0,1,1,0,1,0,1,1,0,1,0,1,0,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,0,1,1,1,0,1,0,1,1,1,0,1,1,0,0,1,1,0,1,1,1,1,0,1,0,1,0,1,0,1,1,1,1,1,0,1,1,0,1,1,0,1,0,0,1,0,1,1,1,0,1,0,1,1,0,0,0,0,0,1,0,1,1,0,1,0,1,0,0,1,0,1,1,0,0,0,0]);
tv1_doubled = np.vstack((tv1_tp4, tv1_tp4))
tv1_parity = np.array([0,1,0,1,0,0,0,1])
eye_8 = np.eye(8, dtype = np.int32)
eye_60 = np.eye(60, dtype = np.int32)

#This is the generator matrix per https://www.ieee802.org/3/df/public/22_10/22_1005/bliss_3df_01_220929.pdf
G1 = np.vstack((eye_60, G.transpose()))
H1 = np.hstack((G.transpose(), eye_8))

pathToMatrix = projectDir + '/bliss_3df_01_220929.mat'
workspace = scipy.io.loadmat(pathToMatrix)

P = workspace['p']

I_120_120 = np.identity(120)
G1 = np.vstack((P,I_120_120))

H1 = np.hstack((np.identity(8), P))

# Check H is a parity matrix for G:
assert(np.all(H1.dot(G1) % 2 ==0))