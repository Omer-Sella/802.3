# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 20:02:48 2024

@author: Omer
"""
from itertools import combinations, combinations_with_replacements
import heapq

def chase2Decoder(hardDecisionDecoder, receivedVector, scores, numberOfLeastProbable):
    """
    Chase 2 decoder takes a hard decision decoder hardDecisionDecoder(receivedVector),
    and engineers a soft input decoder out of it by generating multiple testVectors and
    attempting to decode them. 
    The testVectors are bit flips of the least probable (lowest score) received bits
        
    """
    
    leastProbableIndices = heapq.nsmallest(numberOfLeastProbable, scores, key=None)
    # Should be equivalent to:
    #leastProbableIndices = sorted(scores)[:numberOfLeastProbable]
    # I need to check which is faster.
    
    # We want all test vectors with UP TO numberOfLeastProbable leastProbableIndices, 
    # so we're allowing replacements, but (!), this is not equivalent to FIRST 
    # defining how many we are willing to replace, then find all least probables, 
    # then replace BECAUSE: we are allowing replacement of more probable symbols
    # WITHOUT replacing the lesser ones. It's an implementation choice that should 
    # be measured.
    testIndices = list(combinations_with_replacements(leastProbableIndices, numberOfLeastProbable))
    # Need something like zip
    testVectors = [ x[y] = 1 - x[y] for (x,y) in zip(receivedVector, testIndices)]
    # Omer: how do we deal with decoder failure ?
    decodedVectors = list(map(hardDecisionDecoder, testVectors))
    decodedVector = np.apply_along_axis(np.bincount ...
    return decodedVectors
    