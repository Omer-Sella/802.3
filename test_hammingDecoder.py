# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 17:14:04 2024

@author: Omer
"""
from ieeeConstants import parityMatrix_177_5, g_177_1, IEEE_8023_INT_DATA_TYPE
from ieee8023dj_d0p1 import encode_177_5
from hammingBasics import simpleHammingDecoder, hammingWrapper128_68_128
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
    bestVector, decoderFailure, results = chase2Decoder(wiredDecoder, receivedVector = zeroCodeword, scores = scores, numberOfLeastProbable = 2)
    for r in results:
        print(r)
        
def test_chase2ZeroDataOneOrTwoErrorsLabConditions():
    """
    The purpose of this test is to check that, if just one or two errors occurs AND they are the lowest scoring bits, 
    with all others at 100% confidence, chase2 decodes only that single error 
    """
    def wiredDecoder(x):
        return hammingWrapper128_68_128(x)
    #binaryData = np.random.randint(0, 2, dtype = IEEE_8023_INT_DATA_TYPE)
    #encodedBinaryData
    zeroCodeword = np.zeros(128, dtype = IEEE_8023_INT_DATA_TYPE)
    scores = np.ones(128, dtype = IEEE_8023_INT_DATA_TYPE)
    #scores = np.random.uniform(low = 0, high = 1, size = 128) #
    bestVector, decoderFailure, results = chase2Decoder(wiredDecoder, receivedVector = zeroCodeword, scores = scores, numberOfLeastProbable = 2)
    for r in results:
        print(r)
        
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
    bestVector, decoderFailure, results = chase2Decoder(wiredDecoder, receivedVector = zeroCodeword, scores = scores, numberOfLeastProbable = 2)
    for r in results:
        print(r)
    
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
        bitsIn[indices] = 0
    

    


    