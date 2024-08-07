# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 17:14:04 2024

@author: Omer
"""
from ieeeConstants import parityMatrix_177_5, g_177_1, IEEE_8023_INT_DATA_TYPE
from ieee8023dj_d0p1 import encode_177_5
from hammingBasics import simpleHammingDecoder, hammingWrapper128_68_128, indexToSyndrome, syndromeToIndex
from chase2 import chase2Decoder
import numpy as np
from itertools import combinations

def test_clause_177_1_encoderToParityMatrix_177_5(someBigNumber = 10000):
    def smallEncoder(x):
        return encode_177_5(g_177_1, x)
    #Note that the parity matrix is not in the specification, but we can use Will Bliss's presentation to create it, which we do in constants.
    data = np.random.randint(0,2,(someBigNumber,120), IEEE_8023_INT_DATA_TYPE)
    codewords = list(map(smallEncoder, [data[i,:] for i in range(data.shape[0])]))
    syndromes = [parityMatrix_177_5.dot(codewords[i].transpose()) % 2 for i in range(len(codewords))]
    for s in syndromes:
        assert(np.all(s == 0))
        
def test_detectOneError(someBigNumber = 1000):
    def smallEncoder(x):
        return encode_177_5(g_177_1, x)
    #Note that the parity matrix is not in the specification, but we can use Will Bliss's presentation to create it, which we do in constants.
    data = np.random.randint(0,2,(someBigNumber,120), IEEE_8023_INT_DATA_TYPE)
    codewords = [smallEncoder(data[i,:]) for i in range(someBigNumber)]
    for c in codewords:
        index = np.random.randint(0,68)
        c[index] = 1 - c[index]
    syndromes = [parityMatrix_177_5.dot(codewords[i].transpose()) % 2 for i in range(len(list(codewords)))]
    for s in syndromes:
        assert len(s) == 8
        assert(np.sum(s) != 0)
        
def test_correctOneErrorZeroDataExhaustiveCoverage():
    def smallDecoder(r):
        return simpleHammingDecoder(np.array(parityMatrix_177_5), r) # the constant was given as an np.matrix, and np.array is better.
    encodedZeroData = np.zeros(68, IEEE_8023_INT_DATA_TYPE)
    for index in range(68):
        encodedZeroData[index] = 1
        correctedVector, correctionVector, decoderFailure, syndrome = smallDecoder(encodedZeroData)
        assert (np.sum(correctionVector) == 1)
        assert (correctionVector[index] == 1)
        encodedZeroData[index] = 0
    
    
def test_correctOneErrorRandomDataSimpleHammingDecoder(someBigNumber = 1000):
    def smallEncoder(x):
        return encode_177_5(g_177_1, x)
    def smallDecoder(r):
        return simpleHammingDecoder(np.array(parityMatrix_177_5), r) # the constant was given as an np.matrix, and np.array is better.
    #Note that the parity matrix is not in the specification, but we can use Will Bliss's presentation to create it, which we do in constants.
    data = np.random.randint(0,2,(someBigNumber,120), dtype =  IEEE_8023_INT_DATA_TYPE)
    codewords = [smallEncoder(data[i,:]) for i in range(someBigNumber)]
    indecies = []
    index = np.random.randint(0,68)
    for c in codewords:
        c[index] = 1 - c[index]
        correctedVector, correctionVector, decoderFailure, syndrome = smallDecoder(c)
        assert (np.sum(correctionVector) == 1)
        assert (correctionVector[index] == 1)
        
        
def test_simpleHammingDecoder():
    for i in range(68):
        error = np.zeros(68, dtype = IEEE_8023_INT_DATA_TYPE)
        error[i] = 1
        correctedMessage, correctionVector, decoderFailure, syndrome = simpleHammingDecoder(np.array(parityMatrix_177_5), error)
        assert(np.all(correctedMessage == 0))
        assert(np.all(error == correctionVector))
    

def test_simpleHammingDecoderZeroCodeword():
    error = np.zeros(68, dtype = IEEE_8023_INT_DATA_TYPE)
    correctedMessage, correctionVector, decoderFailure, syndrome = simpleHammingDecoder(np.array(parityMatrix_177_5), error)
    assert(decoderFailure == False)
    assert(np.all(correctedMessage == 0))
    assert(np.all(correctionVector == 0))
    assert(np.all(syndrome == 0))
    
    
def test_chase2ZeroDataNoErrorsFullScore():
    """
    The purpose of this test is to check that, if the all zero codeword is received, 
    with 100% confidence, chase2 agrees that this is the best codeword.
    """
    def wiredDecoder(x):
        return hammingWrapper128_68_128(x)
    zeroCodeword = np.zeros(128, dtype = IEEE_8023_INT_DATA_TYPE)
    scores = np.ones(128, dtype = IEEE_8023_INT_DATA_TYPE)
    bestVector, score, decoderFailure, results = chase2Decoder(wiredDecoder, receivedVector = zeroCodeword, scores = scores, numberOfLeastProbable = 2)
    assert (np.all(bestVector == 0))
    assert (decoderFailure == False)
    for r in results:
        assert (np.all(r[0] == 0))
        assert (r[2] == False)
        
def test_chase2ZeroDataOneErrorLabConditions():
    """
    The purpose of this test is to check that, if just one error occurs AND it is the lowest scoring bits, 
    with all others at 100% confidence, chase2 decodes only that single error 
    """
    def wiredDecoder(x):
        return hammingWrapper128_68_128(x)
    #binaryData = np.random.randint(0, 2, dtype = IEEE_8023_INT_DATA_TYPE)
    #encodedBinaryData
    zeroCodeword = np.zeros(128, dtype = IEEE_8023_INT_DATA_TYPE)
    scores = np.ones(128, dtype = IEEE_8023_INT_DATA_TYPE)
    for i in range(128):
        # Flip one bit and set score to 0
        zeroCodeword[i] = 1 - zeroCodeword[i]
        scores[i] = 0
        bestVector, score, decoderFailure, results = chase2Decoder(wiredDecoder, receivedVector = zeroCodeword, scores = scores, numberOfLeastProbable = 2)    
        assert (np.all(bestVector == 0))
        assert (decoderFailure == False)
        # Reset for next testing
        zeroCodeword[i] = 1 - zeroCodeword[i]
        scores[i] = 1
def test_chase2ZeroDataTwoErrorsLabConditions():
    """
    The purpose of this test is to check that, if just two errors occured AND it they are the lowest scoring bits, 
    with all others at 100% confidence, chase2 decodes only those two errors
    """
    def wiredDecoder(x):
        return hammingWrapper128_68_128(x)
    #binaryData = np.random.randint(0, 2, dtype = IEEE_8023_INT_DATA_TYPE)
    #encodedBinaryData
    zeroCodeword = np.zeros(128, dtype = IEEE_8023_INT_DATA_TYPE)
    scores = np.ones(128, dtype = IEEE_8023_INT_DATA_TYPE)
    for i in list(combinations(range(128),2)):
        # Flip one bit and set score to 0
        zeroCodeword[list(i)] = 1 - zeroCodeword[list(i)]
        scores[list(i)] = 0
        bestVector, score, decoderFailure, results = chase2Decoder(wiredDecoder, receivedVector = zeroCodeword, scores = scores, numberOfLeastProbable = 2)    
        assert (np.all(bestVector == 0))
        assert (decoderFailure == False)
        # Reset for next testing
        zeroCodeword[list(i)] = 1 - zeroCodeword[list(i)]
        scores[list(i)] = 1
        
def test_chase2ZeroDataTwoErrorsFavorableConditions():
    """
    The purpose of this test is to check that if the all zero word is used and there are two bit flips - with lowest scores - 
    chase2  will decode these lowest scoring bits.
    Other bits will have some score between 0 and 1, but higher than the two flipped bits.
    """
    def wiredDecoder(x):
        return hammingWrapper128_68_128(x)
    #binaryData = np.random.randint(0, 2, dtype = IEEE_8023_INT_DATA_TYPE)
    #encodedBinaryData
    zeroCodeword = np.zeros(128, dtype = IEEE_8023_INT_DATA_TYPE)
    scores = np.ones(128, dtype = IEEE_8023_INT_DATA_TYPE)
    #scores = np.random.uniform(low = 0, high = 1, size = 128) #
    for i in list(combinations(range(128),2)):
        # Flip one bit and set score to 0
        zeroCodeword[list(i)] = 1 - zeroCodeword[list(i)]
        scores[list(i)] = np.random.uniform(low = 0, high = 0.8, size = 2)
        bestVector, score, decoderFailure, results = chase2Decoder(wiredDecoder, receivedVector = zeroCodeword, scores = scores, numberOfLeastProbable = 2)    
        assert (np.all(bestVector == 0))
        assert (decoderFailure == False)
        # Reset for next testing
        zeroCodeword[list(i)] = 1 - zeroCodeword[list(i)]
        scores[list(i)] = 1
    
def test_hammingReportFailure():
    # I need to fix this test, basically testing patterns overwhich hamming should declare failure.
    def wiredDecoder(x):
        return hammingWrapper128_68_128(x)
    pass
    #bitsIn = np.zeros(128, IEEE_8023_INT_DATA_TYPE)
    #pairs = combinations(range(128), 2)
    #for p in pairs:
    #    indices = list(p)
    #    print(indices)
    #    bitsIn[indices] = 1
    #    print(bitsIn)
    #    correctedMessage128, correctionVector128, decoderFailure, syndromes = wiredDecoder(bitsIn)
        #assert(decoderFailure == True)
    #        bitsIn[indices] = 0
    

def test_bliss_3df_01_220929():
    #I_120_120 = np.identity(120)
    #G1 = np.vstack((P,I_120_120))
    
    #H1 = np.hstack((np.identity(8), P))
    
    # Check H is a parity matrix for G:
    #assert(np.all(H1.dot(G1) % 2 ==0))
    import scipy
    import os
    import sys
    projectDir = os.environ.get('IEEE8032DJ')
    #You don't have to define an environment variable, but then you have to give a path to the project here
    if projectDir == None: 
         projectDir = "c:/users/omer/802.3/"
    sys.path.insert(1, projectDir)
    pathToMatrix = projectDir + '//bliss_3df_01_220929.mat'
    matrices = scipy.io.loadmat(projectDir + '/bliss_3df_01_220929.mat')
    p = matrices['p']
    parityMatrixFrom_bliss_3df_01_220929 = matrices['h']
    
    
    I_120_120 = np.identity(120)
    I_8_8 = np.identity(8)
    
    G = np.hstack(( p.transpose(), I_120_120))
    
    # This is the parity matrix presented as Canonic in the presentation
    h1 = np.vstack((I_8_8, p.transpose()))
    h1 = h1.transpose()
    
    sanityCheck = G.dot(parityMatrixFrom_bliss_3df_01_220929.transpose()) % 2
    assert(np.all(sanityCheck == 0))

def test_indexToSyndromeToIndex():
    for i in np.arange(0,128,1):
        string = np.binary_repr(i, 8)
        array = [int(s) for s in string]
        assert(np.all(np.flipud(syndromeToIndex(indexToSyndrome(array)))[0:7] == array[1:8]))




def test_hammingWrapper128_68_128ZeroCodeword():
    data = np.zeros(128, dtype = IEEE_8023_INT_DATA_TYPE)
    correctedMessage128, correctionVector128, decoderFailure, syndromes = hammingWrapper128_68_128(data)
    assert (np.all(correctedMessage128 == 0))
    assert (np.all(correctedMessage128 == data)) # redundant but still ...
    assert (np.all(correctionVector128 == 0))
    assert (decoderFailure == False)
    assert (np.all(syndromes == 0))
    
