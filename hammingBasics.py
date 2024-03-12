# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 19:09:31 2023

@author: Megatron
"""

import scipy.io
import numpy as np
from ieeeConstants import *
BLISS_LOGICAL_DATA_TYPE = np.int32
#please refer to https://www.ieee802.org/3/df/public/22_10/22_1005/bliss_3df_01_220929.pdf
INDEX_TO_NUMBER = np.array([ 0, 1, 2, 4, 8, 16, 32, 64]).transpose()
projectDir = 'd:/802.3/'
matrices = scipy.io.loadmat(projectDir + 'bliss_3df_01_220929.mat')
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

def bliss_3df_01_220929_syndrom_decoder(slicedReceivedMessage, H):

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
    syndrome = np.squeeze(np.asarray(slicedReceivedMessage.dot(H)))
    correctionVector = np.zeros(slicedReceivedMessage.shape[0], dtype = IEEE_8023_INT_DATA_TYPE)
    index = 0
    decoderFailure = False
    if np.all(syndrome == 0):
        pass
    else:
        found = False
        while (index < H.shape[0]) and (not found):
            if np.all(H[index , :] == syndrome):
                found = True
            else:
                index = index + 1
            
        if index >= H.shape[0]:
                decoderFailure = True
        else:
                correctionVector[index] = 1
        
    correctedMessage = (slicedReceivedMessage + correctionVector) %2
    return  correctedMessage, correctionVector, decoderFailure
                
def test_simpleHammingDecoder():
    for i in range(68):
        error = np.zeros(68, dtype = IEEE_8023_INT_DATA_TYPE)
        error[i] = 1
        correctedMessage, correctionVector, decoderFailure = simpleHammingDecoder(H1.transpose(), error)
        assert(np.all(correctedMessage == 0))
        assert(np.all(error == correctionVector))
    return 'OK'
