# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 19:09:31 2023

@author: Megatron
"""

import scipy.io
import numpy as np
from ieeeConstants import *
from ieeeConstants import parityMatrix_177_5
BLISS_LOGICAL_DATA_TYPE = np.int32
#please refer to https://www.ieee802.org/3/df/public/22_10/22_1005/bliss_3df_01_220929.pdf
INDEX_TO_NUMBER = np.array([ 0, 1, 2, 4, 8, 16, 32, 64]).transpose()
import os
import sys
projectDir = os.environ.get('IEEE8032DJ')
#You don't have to define an environment variable, but then you have to give a path to the project here
if projectDir == None: 
     projectDir = "c:/users/omer/802.3/"
sys.path.insert(1, projectDir)
pathToMatrix = projectDir + '//bliss_3df_01_220929.mat'
#workspace = scipy.io.loadmat(pathToMatrix)


#P = workspace['p']

#I_120_120 = np.identity(120)
#G1 = np.vstack((P,I_120_120))

#H1 = np.hstack((np.identity(8), P))

# Check H is a parity matrix for G:
#assert(np.all(H1.dot(G1) % 2 ==0))
matrices = scipy.io.loadmat(projectDir + '/bliss_3df_01_220929.mat')
p = matrices['p']
h = matrices['h']


I_120_120 = np.identity(120)
I_8_8 = np.identity(8)

G = np.hstack(( p.transpose(), I_120_120))

# This is the parity matrix presented as Canonic in the presentation
h1 = np.vstack((I_8_8, p.transpose()))
h1 = h1.transpose()


sanityCheck = G.dot(h.transpose()) % 2


#if (np.all(sanityCheck == 0)):
#    print('Seems that h is a parity for G.')
#else:
#    print('Check your parity - generator pair')
    

#Observation - to get a coded message, multiply a message on the left by G on the right: codedMessage = message.dot(G)

def bliss_3df_01_220929_syndrom_decoder(H, slicedReceivedMessage):

    #S = H.dot(slicedReceivedMessage) %2
    S = slicedReceivedMessage.dot(H.transpose()) %2
    print(S)
    index = np.zeros(8, dtype = BLISS_LOGICAL_DATA_TYPE)
    index[0] = S[0]
    index[1] = S[1] 
    index[2] = S[0] * S[1] + S[2]  %2
    index[3] = (1 - S[0]) * (1 - S[1]) * S[2] + S[3] %2
    index[4] = S[0] * (1 - S[1]) * S[2]  + S[4] %2
    index[5] = (1 - S[0]) * S[1] * S[2] + S[5] %2
    index[6] = S[0] * S[1] * (1 - S[2]) + S[6] %2
    
    integerIndex = index.dot(INDEX_TO_NUMBER)
    
    correctionVector = np.zeros(slicedReceivedMessage.shape[0], dtype = BLISS_LOGICAL_DATA_TYPE)
    correctionVector[integerIndex] = 1
    
    #Create a corrected message
    correctedMessage = slicedReceivedMessage + correctionVector %2
    
    #Check if the corrected message is a codeword, if not, report decoder failure.
    decoderFailure = np.all(H.dot(correctedMessage) == 0)
    
    return correctedMessage, correctionVector, decoderFailure

def simpleHammingDecoder(H, slicedReceivedMessage):
    """
    H is parity matrix, assumed to have unique binary columns, i.e.: it can correct one bit flip.
    slicedReceivedMessage is assumed to be a binary vector of length compatible to the dimensions of H 
    
    decoderFailue - is a flag that indicates whether the decoder found a correction (or if non was needed)
    correctionVector - should be a binary vector of at most one nonzero
    correctedMessage - a binary vector after correction was applied, i.e.: (slicedReceivedMessage + correctionVector) %2
    
    """
    correctionVector = np.zeros(len(slicedReceivedMessage), dtype = IEEE_8023_INT_DATA_TYPE)
    decoderFailure = False
    #syndrome = slicedReceivedMessage.dot(H)
    syndrome = H.dot(slicedReceivedMessage.transpose()) % 2
    if np.all(syndrome == 0):
        decoderFailure = False
    else:
    #syndrome = np.squeeze(np.asarray(slicedReceivedMessage.dot(H)))
        correctionVector = [np.all(syndrome == H[:,i]) for i in range(H.shape[1])]
        if np.all(correctionVector == 0):
            decoderFailure = True
    correctedMessage = (slicedReceivedMessage + correctionVector) %2
    return  correctedMessage, correctionVector, decoderFailure

def hammingWrapper(bitsIn):
    """
    This function XORs bit pairs from 0 to 120, and leaves rightmost 8 bits untouched, then applies Hamming
    hammingDecoder68BitFunction is assumed to be a Hamming decoder, already set to some parity matrix H
    """
    def hammingDecoder68BitFunction(x):
        return simpleHammingDecoder(parityMatrix_177_5, x)
    xored = np.zeros(68, dtype = IEEE_8023_INT_DATA_TYPE)
    xored[0:60] =  [(bitsIn[2*k] + bitsIn[2 * k + 1]) %2 for k in range(60) ]
    xored[60:68] =  bitsIn[120:128]
    correctedMessage, correctionVector, decoderFailure = hammingDecoder68BitFunction(xored)
    return correctedMessage, correctionVector, decoderFailure
    

def hamming_177_1_using_polynomials():
    #The key here was to find a polynomial that generates table 177_1, and I think it is  x^8 + x^7 + x^6 + x^5 + x^2 + x^1 + 1, except the notation is right to left, i.e.: the leading coefficient is on the right, and the free element (0 or 1 ) is on the left
    pass
                

