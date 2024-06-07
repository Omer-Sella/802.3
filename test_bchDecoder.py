# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 14:18:16 2024

@author: Omer
"""
import os, sys
projectDir = os.environ.get('IEEE8032DJ')
#You don't have to define an environment variable, but then you have to give a path to the project here
if projectDir == None: 
     projectDir = "c:/Users/Omer/802.3/"
sys.path.insert(0, projectDir)
reedSolomonProjectDir = os.environ.get('REEDSOLOMON')
if reedSolomonProjectDir == None: 
     reedSolomonProjectDir = "c:/users/omer/reedSolomon/reedSolomon/"
sys.path.insert(0, reedSolomonProjectDir)
from arithmetic import binaryFieldElement as galoisElement
from arithmetic import polynomial as polynomialClass
from arithmetic import generateExponentAndLogTables, gf128, generateExponentAndLogTables
from ieee8023dj_d0p1 import bchEncoder
from bchDecoder import bchDecoder
import numpy as np


def test_bchDecoder():
    zeroData = np.zeros(110)
    encodedZeroData = bchEncoder(zeroData)
    eD, _ =  generateExponentAndLogTables()
    #print("Done generating log and exp dictionaries.")
    correctedVector, correctionVector, errorLocatorX = bchDecoder( receivedBinaryVecotor = encodedZeroData, exponentDictionary = eD, numberOfPowers = 16, codewordLengthMaximal = 127)
    assert (np.all(correctedVector == 0))
    
def test_bchDecoder_single_bit_flip():
    zeroData = np.zeros(110)
    eD, _ =  generateExponentAndLogTables()
    encodedZeroData = bchEncoder(zeroData)
    for i in range(len(encodedZeroData)):
        encodedZeroData[i] = 1
        correctedVector, correctionVector, errorLocatorX = bchDecoder( receivedBinaryVecotor = encodedZeroData, exponentDictionary = eD, numberOfPowers = 16, codewordLengthMaximal = 127)
        print("Error locator polynomial: ")
        errorLocatorX.printValues()
        print("corrected vector:")
        print(correctedVector)
        print("correction vector:")
        print(correctionVector)
        print("****")
        assert correctionVector[i] == 1
        encodedZeroData[i] = 0

def test_bchDecoder_single_bit_flip_order_check():
    zeroData = np.zeros(110)
    eD, _ =  generateExponentAndLogTables()
    encodedZeroData = bchEncoder(zeroData)
    for i in range(len(encodedZeroData)):
        encodedZeroData[i] = 1
        correctedVector, correctionVector, errorLocatorX = bchDecoder( receivedBinaryVecotor = encodedZeroData, exponentDictionary = eD, numberOfPowers = 16, codewordLengthMaximal = 127)
        errorLocatorX.printValues()
        assert errorLocatorX.order() == 1
        encodedZeroData[i] = 0
        
def coverage_retrace_bug_error_in_first_coordinate(index):
    zeroData = np.zeros(110)
    eD, _ =  generateExponentAndLogTables()
    encodedZeroData = bchEncoder(zeroData)
    encodedZeroData[index] = 1
    # Notice that the decoder needs to produce the error locator polynomial eX for this coverage !
    correctedVector, correctionVector, eX = bchDecoder( receivedBinaryVecotor = encodedZeroData, exponentDictionary = eD, numberOfPowers = 16, codewordLengthMaximal = 127)
    for i in eD.keys():
        print(eX.at(gf128(eD[i])).getValue())
    return correctedVector, correctionVector, eX
    

def coverage_example_8_8_tkmoon():
    pass

def test_connection_polynomial_for_two_errors_explicit_calculation():   
    zeroData = np.zeros(110)
    eD, _ =  generateExponentAndLogTables()
    encodedZeroData = bchEncoder(zeroData)
    encodedZeroData[1] = 1
    encodedZeroData[10] = 1
    # Notice that the decoder needs to produce the error locator polynomial eX for this coverage !
    correctedVector, correctionVector, eX = bchDecoder( receivedBinaryVecotor = encodedZeroData, exponentDictionary = eD, numberOfPowers = 16, codewordLengthMaximal = 127)
    lambda0 = gf128(1)
    lambda1 = syndromeAsGf128[0]
    #lambda2 = (syndromeAsGf128[3] + (syndromeAsGf128[1] * syndromeAsGf128[1] * syndromeAsGf128[1]) / syndromeAsGf128[0])
    explicitConnectionX = polynomial(coefficients = [lambda1, lambda0])
    assert explicitConnectionX == cX