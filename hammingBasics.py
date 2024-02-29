# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 19:09:31 2023

@author: Megatron
"""

import scipy.io
import numpy as np
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


    
def test_isSingleErrorCorrecting(H):
    oneHot = np.zeros(H.shape[1])
    syndromeList = []
    for i in range(H.shape[1]):
        oneHot[i] = 1
        syndromes = h.dot(oneHot) %2
        syndromeList.append(syndromes)
        oneHot[i] = 0
    assert (len(np.unique(syndromeList, axis = 0)) == H.shape[1])
    return syndromeList

# def test_hardDecoder(decoderFunction = simpleHammingDecoder):
#     message = np.random.randint(0,2,120)
#     codedMessage = message.dot(G) %2

#     #Observation - to check / get syndromes, multiply the parity on the right by a coded message on the left: syndromes = h.dot(codedMessage)

#     syndromes = h1.dot(codedMessage) %2

#     if (np.all(syndromes) == 0):
#         print('Valid codeword received, no decoding needed.')
#     else:
#         print('InValid vector. Decoding needed.')
    
#     for i in range(len(codedMessage)):
#         print(i)
#         codedMessage[i] = 1 - codedMessage[i]
#         cm,cv, df = decoderFunction(codedMessage, h)
#         if ( (cv[i] == 1) and (np.sum(cv) == 1) and (df == False)):
#             print('OK for a single error at index ' + str(i))
#         else:
#             print('Single error at index ' + str(i) + 'issue.')
    
#     return 'OK'

