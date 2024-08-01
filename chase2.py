# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 20:02:48 2024

@author: Omer
"""
from itertools import combinations, combinations_with_replacement
import heapq
import copy
import numpy as np
from ieeeConstants import IEEE_8023_INT_DATA_TYPE

def generateTestVectors(template, indices, replacements):
    testVector = copy.copy(template)
    testVector[indices] = replacements[indices]
    return testVector

def chase2Decoder(hardDecisionDecoder, receivedVector, scores, numberOfLeastProbable = 2, scoreDecoderCorrection = False):
    """
    Chase 2 decoder takes a hard decision decoder hardDecisionDecoder(receivedVector),
    and engineers a soft input decoder out of it by generating multiple testVectors and
    attempting to decode them. 
    The testVectors are bit flips of the least probable (lowest score) received bits
    
    Assumptions:
        1. We work with bits, not PAM4 symbols (or anything else), therefore:
        2. A replacement / swap is simply a bit flip, and:
        3. If bit b has score s then bit 1-b has score 1-s
        
        The algorithm is best described by Cathy Liu in https://www.signalintegrityjournal.com/articles/3405-200-gbps-ethernet-forward-error-correction-fec-analysis        
        1. Form the HD received sequence z from r and assign a reliability value to each symbol of z
        2. Generate the error patterns in E one at a time, possible in likelihood order. For each error pattern e in E, form the test patterns z+e
        3. Decode each test pattern into a codeword using HD decoder
        4. Compute the soft decision decoding metric for each generated candidate codeword
        5. Select the candidate codeword with the best metric as the coded solution.
        
        Where the implementation of 5 in this work is as follows:
        The score for a vector is the sum of its bit scores (decoder failure is reflected by setting the score to 0)
        
    """
    
    
    #leastProbable = heapq.nsmallest(numberOfLeastProbable, scores, key=None)
    #leastProbable = sorted(scores)[:numberOfLeastProbable]
    leastProbable = heapq.nsmallest(numberOfLeastProbable, scores)
    indices = np.array([], dtype = np.int32)
    for l in leastProbable:
        indices = np.hstack((indices , np.int32(np.squeeze(np.where(scores == l)))))
    indices = np.unique(indices)
    #print(leastProbable)
    #print(list(indices))
    
       
    # We want all test vectors with UP TO numberOfLeastProbable leastProbableIndices, 
    # so we're allowing replacements, but (!), this is not equivalent to FIRST 
    # defining how many we are willing to replace, then find all least probables, 
    # then replace BECAUSE: we are allowing replacement of more probable symbols
    # WITHOUT replacing the lesser ones. It's an implementation choice that should 
    # be measured
    testIndices = combinations_with_replacement(indices, numberOfLeastProbable)
    #score = 0
    bestScore = 0
    bestVector = np.array([])
    decoderFailure = True
    results = []
    score = np.sum(scores)
    initialScore = np.sum(scores)
    for t in testIndices:
        tAsList = list(t)
        newTestVector = copy.deepcopy(receivedVector)
        newTestVector[tAsList] = 1 - receivedVector[tAsList]
        testCorrectedMessage, testCorrectionVector, testDecoderFailure, testSyndromes = hardDecisionDecoder(newTestVector)
        if testDecoderFailure == False:
            # At least one test vector managed to get through
            decoderFailure = False
            # Update the score of the current test vector. 
            # The score is made up of the sum of scores, 
            # minus the scores of the bits we flipped (either for test or by the deocder, with possible overlap) 
            # plus 1-score for those bits that we flipped (again, as test or flipped by the decoder)
            
            if scoreDecoderCorrection:
                score = initialScore - np.sum(scores[tAsList]) + np.sum(1 - scores[tAsList]) - np.sum(scores[testCorrectionVector]) + np.sum(1 - scores[testCorrectionVector])
            else:
                score = initialScore - np.sum(scores[tAsList]) + np.sum(1 - scores[tAsList])
            if score > bestScore:
                bestVector = testCorrectedMessage
                bestScore = score
            results.append([testCorrectedMessage, score, testDecoderFailure])
    
    return bestVector, score, decoderFailure, results

