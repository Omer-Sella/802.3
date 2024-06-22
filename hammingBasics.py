# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 19:09:31 2023

@author: Megatron
"""

import scipy.io
import numpy as np
from ieeeConstants import *
from ieeeConstants import parityMatrix_177_5, IEEE_8023_INT_DATA_TYPE
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

# #I_120_120 = np.identity(120)
# #G1 = np.vstack((P,I_120_120))

# #H1 = np.hstack((np.identity(8), P))

# # Check H is a parity matrix for G:
# #assert(np.all(H1.dot(G1) % 2 ==0))
# matrices = scipy.io.loadmat(projectDir + '/bliss_3df_01_220929.mat')
# p = matrices['p']
# parityMatrixFrom_bliss_3df_01_220929 = matrices['h']


# I_120_120 = np.identity(120)
# I_8_8 = np.identity(8)

# G = np.hstack(( p.transpose(), I_120_120))

# # This is the parity matrix presented as Canonic in the presentation
# h1 = np.vstack((I_8_8, p.transpose()))
# h1 = h1.transpose()


# sanityCheck = G.dot(parityMatrixFrom_bliss_3df_01_220929.transpose()) % 2


#if (np.all(sanityCheck == 0)):
#    print('Seems that h is a parity for G.')
#else:
#    print('Check your parity - generator pair')
    

#Observation - to get a coded message, multiply a message on the left by G on the right: codedMessage = message.dot(G)

def bliss_3df_01_220929_syndrom_decoder(H, slicedReceivedMessage):

    #S = H.dot(slicedReceivedMessage) %2
    S = H.dot(slicedReceivedMessage) %2
    S = S[::-1] # Bliss defines S8 to be the top
    
    print(S)
    index = np.zeros(8, dtype = IEEE_8023_INT_DATA_TYPE)
    index[1] = S[1-1]
    index[2] = S[2-1] 
    index[3] = ((S[1-1] * S[2-1]) + S[3-1])  %2
    index[4] = (((1 - S[1-1]) * (1 - S[2-1]) * S[3-1]) + S[4-1]) %2
    index[5] = (S[1-1] * (1 - S[2-1]) * S[3-1]  + S[5-1]) %2
    index[6] = ((1 - S[1]) * S[2-1] * S[3-1] + S[6-1]) %2
    index[7] = S[2-1] * S[2-1] * (((1 - S[3-1]) + S[7-1]) %2)
    
    print(index)
    integerIndex = index.dot(INDEX_TO_NUMBER)
    
    correctionVector = np.zeros(slicedReceivedMessage.shape[0], dtype = BLISS_LOGICAL_DATA_TYPE)
    correctionVector[integerIndex] = 1
    
    #Create a corrected message
    correctedMessage = (slicedReceivedMessage + correctionVector) %2
    
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
    decoderFailure = True
    #syndrome = slicedReceivedMessage.dot(H)
    syndrome = H.dot(slicedReceivedMessage.transpose()) % 2
    if np.all(syndrome == 0):
        decoderFailure = False
    else:
    #syndrome = np.squeeze(np.asarray(slicedReceivedMessage.dot(H)))
        correctionVector = [np.all(syndrome == H[:,i]) for i in range(H.shape[1])]
        if np.all(correctionVector == 0):
            decoderFailure = False
    correctedMessage = (slicedReceivedMessage + correctionVector) %2
    return  correctedMessage, correctionVector, decoderFailure, syndrome

def hammingWrapper128_68_128(bitsIn):
    """
    This function XORs bit pairs from 0 to 120, and leaves rightmost 8 bits untouched, then applies Hamming
    hammingDecoder68BitFunction is assumed to be a Hamming decoder, already set to some parity matrix H
    
    If the decoding succeeds, the result has to be run again through the decoder, since we need to figure out which of the bits from the xored bit-pairs was flipped
    """
    def hammingDecoder68BitFunction(x):
        return simpleHammingDecoder(parityMatrix_177_5, x)
    xored = np.zeros(68, dtype = IEEE_8023_INT_DATA_TYPE)
    xored[0:60] =  [(bitsIn[2*k] + bitsIn[2 * k + 1]) %2 for k in range(60) ]
    xored[60:68] =  bitsIn[120:128]
    correctedMessage128 = np.ones(128, IEEE_8023_INT_DATA_TYPE)
    correctionVector128 = np.zeros(128, IEEE_8023_INT_DATA_TYPE)
    correctedMessage, correctionVector, decoderFailure, syndromes = hammingDecoder68BitFunction(xored)
    print(xored)
    #print(correctedMessage)
    #print(correctionVector)
    print(decoderFailure)
    print(syndromes)
    # Now we need to reverse engineer which of the original bits the hamming-corrected-bit corresponds to.
    # We do this in a very simple minded way - try to correct both, and see which one comes back as a valid codeword:
    if not decoderFailure and not (np.all(syndromes == 0)): #I.e.: Hamming thinks there is (up to) one error and it found it, meaning the correctionVector is one hot or all zero
        print(correctionVector)    
        assert(np.sum(correctionVector) == 1) # Safety - remove when you think everything works. I need to add a test for this.
        oneHotIndex = np.where(correctionVector == 1)[0]
        if oneHotIndex >= 60:
            correctionVector128 = np.zeros(128, IEEE_8023_INT_DATA_TYPE) 
            correctionVector128[oneHotIndex] = 1
            correctedMessage128 = bitsIn + correctionVector128
        else:
            # Try flipping location 2 * oneHotIndex and see what comes back
            correctionVector128 = np.zeros(128, IEEE_8023_INT_DATA_TYPE) 
            correctionVector128[2 * oneHotIndex] = 1
            correctedVector128 = bitsIn + correctionVector128
            xored = np.zeros(68, dtype = IEEE_8023_INT_DATA_TYPE)
            xored[0:60] =  [(correctedVector128[2*k] + correctedVector128[2 * k + 1]) %2 for k in range(60) ]
            xored[60:68] =  correctedVector128[120:128]
            correctedMessage0, correctionVector0, decoderFailure0, syndromes0 = hammingDecoder68BitFunction(xored)
            if not decoderFailure0 and np.all(syndromes0 == 0):
                assert(np.all(correctedMessage0 == correctedMessage))
                correctedMessage128 = correctedVector128
            else:
                #It must be the other option, i.e. 2*oneHotIndex + 1
                correctionVector128 = np.zeros(128, IEEE_8023_INT_DATA_TYPE) 
                correctionVector128[2 * oneHotIndex + 1] = 1
                correctedVector128 = bitsIn + correctionVector128
                xored = np.zeros(68, dtype = IEEE_8023_INT_DATA_TYPE)
                xored[0:60] =  [(correctedVector128[2*k] + correctedVector128[2 * k + 1]) %2 for k in range(60) ]
                xored[60:68] =  correctedVector128[120:128]
                correctedMessage1, correctionVector1, decoderFailure1, syndromes1 = hammingDecoder68BitFunction(xored)
                assert(np.all(correctedMessage1 == correctedMessage))
                correctedMessage128 = correctedVector128
    
    return correctedMessage128, correctionVector128, decoderFailure, syndromes
    

def hamming_177_1_using_polynomials():
    #The key here was to find a polynomial that generates table 177_1, and I think it is  x^8 + x^7 + x^6 + x^5 + x^2 + x^1 + 1, except the notation is right to left, i.e.: the leading coefficient is on the right, and the free element (0 or 1 ) is on the left
    pass

def indexToSyndrome(indexArray):
    i = np.flipud(indexArray)
    s = np.zeros(8)
    s[1] = i[1-1]
    s[2] = i[2-1]
    s[3] = (i[3-1] + i[1-1]*i[2-1]) %2
    s[4] = (i[4-1] + ((1+i[1-1])%2)*((1+i[2-1])%2)*(i[3-1] + i[1-1]*i[2-1]))%2
    s[5] = (i[5-1] +  i[1-1]*((1+i[2-1])%2)*((i[3-1]+i[1-1]*i[2-1])%2) )%2
    s[6] = (i[6-1] + ((1+i[1-1])%2)*i[2-1]*( (i[3-1]+i[1-1]*i[2-1])%2))%2
    s[7] = (i[7-1] + i[1-1]*i[2-1]*((1+i[3-1]+i[1-1]*i[2-1])%2))%2
    return s

def syndromeToIndex(s):
    i = np.zeros(8)
    i[1]= s[1] 
    i[2]= s[2]
    i[3]= (s[1]*s[2] + s[3])%2
    i[4]= (((1+s[1])%2)*((1+s[2])%2)*s[3] + s[4])%2
    i[5]= (s[1]*((1+s[2])%2)*s[3] + s[5])%2
    i[6]= ( ((1+s[1])%2)*s[2]*s[3] + s[6])%2
    i[7]= (s[1]*s[2]*((1+s[3])%2) + s[7])%2
    return i

def generateHammingMatrixFromBlissEquations():
    """
        i(1)=S(1) 
        i(2)=S(2)
        i(3)=xor(( S(1) & S(2) ), S(3))           = s(1)*s(2) + s(3)
        i(4)=xor( (~S(1) & ~S(2) & S(3) ), S(4))  = (1+s(1))*(1+s(2))*s(3) + s(4)
        i(5)=xor( (S(1) & ~S(2) & S(3) ), S(5))   = s(1)*(1+s(2))*s(3) + s(5)
        i(6)=xor( (~S(1) & S(2) & S(3) ), S(6))   = (1+s(1))*s(2)*s(3) + s(6)
        i(7)=xor( (S(1) & S(2) & ~S(3) ), S(7))   = s(1)*s(2)*(1+s(3)) + s(7)

        s(1) = i(1)
        s(2) = i(2)
        s(3) = i(3) + s(1)*s(2)                     = i(3) + i(1)*i(2)
        s(4) = i(4) + (1+s(1))*(1+s(2))*s(3)        = i(4) + (1+i(1))*(1+i(2))*(i(3) + i(1)*i(2))
        s(5) = i(5) + s(1)*(1+s(2))*s(3)            = i(5) + i(1)*(1+i(2))*(i(3)+i(1)*i(2))
        s(6) = i(6) + (1+s(1))*s(2)*s(3)            = i(6) + (1+i(1))*i(2)*(i(3)+i(1)*i(2))
        s(7) = i(7) + s(1)*s(2)*(1+s(3))            = i(7) + i(1)*i(2)*(1+i(3)+i(1)*i(2))


    """
    #1. Generate all the numbers from 1 to 128, convert to binary and then to array.
    #2. Then convert the index array to a syndrome, which is a column in the parity matrix
    indices = []
    columns = []
    for i in np.arange(0,128,1):
        string = np.binary_repr(i, 8)
        array = [int(s) for s in string]
        indices.append(array)
        column = indexToSyndrome(array)
        columns.append(column)
    
    return indices, np.array(columns)
    
    