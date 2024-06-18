# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 17:14:04 2024

@author: Omer
"""
from ieeeConstants import parityMatrix_177_5, g_177_1, IEEE_8023_INT_DATA_TYPE
from ieee8023dj_d0p1 import encode_177_5
from hammingBasics import simpleHammingDecoder, hammingWrapper
from chase2 import chase2Decoder
import numpy as np

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
        correctedVector, correctionVector, decoderFailure = smallDecoder(encodedZeroData)
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
        correctedVector, correctionVector, decoderFailure = smallDecoder(c)
        assert (np.sum(correctionVector) == 1)
        assert (correctionVector[index] == 1)
        
        
def test_simpleHammingDecoder():
    for i in range(68):
        error = np.zeros(68, dtype = IEEE_8023_INT_DATA_TYPE)
        error[i] = 1
        correctedMessage, correctionVector, decoderFailure = simpleHammingDecoder(np.array(parityMatrix_177_5), error)
        assert(np.all(correctedMessage == 0))
        assert(np.all(error == correctionVector))
    

def test_simpleHammingDecoderZeroCodeword():
    error = np.zeros(68, dtype = IEEE_8023_INT_DATA_TYPE)
    correctedMessage, correctionVector, decoderFailure = simpleHammingDecoder(np.array(parityMatrix_177_5), error)
    assert(decoderFailure == False)
    assert(np.all(correctedMessage == 0))
    assert(np.all(correctionVector == 0))
    
    
def test_chase2():
    def wiredDecoder(x):
        return hammingWrapper(x)
    #binaryData = np.random.randint(0, 2, dtype = IEEE_8023_INT_DATA_TYPE)
    #encodedBinaryData
    zeroCodeword = np.zeros(128, dtype = IEEE_8023_INT_DATA_TYPE)
    scores = np.ones(128, dtype = IEEE_8023_INT_DATA_TYPE)
    correctedMessage, correctionVector, decoderFailure = hammingWrapper(zeroCodeword)
    
    r = chase2Decoder(hardDecisionDecoder = wiredDecoder, receivedVector = zeroCodeword, scores = scores)
    print(r)
    
    
    


    