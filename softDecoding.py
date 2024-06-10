# -*- coding: utf-8 -*-
"""
Created on Sat Jun  8 21:09:08 2024

@author: Omer
"""
import numpy as np
BINARY_DTYPE = np.int32

def parityToFSM(parityMatrix, statesColumns = True):
    # parityMatrix is assumed to be a binary matrix (having 0s and 1s only).
    # Its columns will be the states
    if statesColumns:
        states = np.arange(0,parityMatrix.shape[1])
        # Safety: check that the columns of the parity matrix are unique:
        assert np.all(np.unique(parityMatrix, axis = 1 ) == parityMatrix)
    else:
        states = np.arange(0,parityMatrix.H1.shape[0])
        # Safety: check that the rows of the parity matrix are unique:
        assert np.all(np.unique(parityMatrix, axis = 0 ) == parityMatrix)
    
    # Build transition matrix
    transitionMatrix = np.zeros((len(states), len(states)), dtype = BINARY_DTYPE)
    #transitionMatrix[i,j] == 1 if 
    return transitionMatrix, emissionMatrix
    