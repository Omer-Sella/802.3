# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 20:25:38 2024

@author: Omer
"""
from ieeeConstants import *

# Obsolete
# def simpleHammingDecoder(H, slicedReceivedMessage):
#     syndrome = np.squeeze(np.asarray(slicedReceivedMessage.dot(H)))
#     correctionVector = np.zeros(slicedReceivedMessage.shape[0], dtype = IEEE_8023_INT_DATA_TYPE)
#     index = 0
#     decoderFailure = False
#     if np.all(syndrome == 0):
#         pass
#     else:
#         found = False
#         while (index < H.shape[0]) and (not found):
#             if np.all(H[index , :] == syndrome):
#                 found = True
#             else:
#                 index = index + 1
            
#         if index >= H.shape[0]:
#                 decoderFailure = True
#         else:
#                 correctionVector[index] = 1
        
#     correctedMessage = (slicedReceivedMessage + correctionVector) %2
#     return  correctedMessage, correctionVector, decoderFailure
                

