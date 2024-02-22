# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 09:44:44 2024

@author: Omer
"""
import numpy as np
G = np.matrix([
    [1,0,0,1,0,1,0,0],
     [0,1,0,0,1,0,1,0],
 [0,0,1,0,0,1,0,1],
 [1,1,0,0,1,0,1,1],
 [1,0,1,1,1,1,0,0],
 [0,1,0,1,1,1,1,0],
 [0,0,1,0,1,1,1,1],
 [1,1,0,0,1,1,1,0],
 [0,1,1,0,0,1,1,1],
 [1,1,1,0,1,0,1,0],
 [0,1,1,1,0,1,0,1],
 [1,1,1,0,0,0,1,1],
 [1,0,1,0,1,0,0,0],
 [0,1,0,1,0,1,0,0],
 [0,0,1,0,1,0,1,0],
 [0,0,0,1,0,1,0,1],
 [1,1,0,1,0,0,1,1],
 [1,0,1,1,0,0,0,0],
 [0,1,0,1,1,0,0,0],
 [0,0,1,0,1,1,0,0],
 [0,0,0,1,0,1,1,0],
 [0,0,0,0,1,0,1,1],
 [1,1,0,1,1,1,0,0],
 [0,1,1,0,1,1,1,0],
 [0,0,1,1,0,1,1,1],
 [1,1,0,0,0,0,1,0],
 [0,1,1,0,0,0,0,1],
 [1,1,1,0,1,0,0,1],
 [1,0,1,0,1,1,0,1],
 [1,0,0,0,1,1,1,1],
 [1,0,0,1,1,1,1,0],
 [0,1,0,0,1,1,1,1],
 [1,1,1,1,1,1,1,0],
 [0,1,1,1,1,1,1,1],
 [1,1,1,0,0,1,1,0],
 [0,1,1,1,0,0,1,1],
 [1,1,1,0,0,0,0,0],
 [0,1,1,1,0,0,0,0],
 [0,0,1,1,1,0,0,0],
 [0,0,0,1,1,1,0,0],
 [0,0,0,0,1,1,1,0],
 [0,0,0,0,0,1,1,1],
[1,1,0,1,1,0,1,0],
[0,1,1,0,1,1,0,1],
[1,1,1,0,1,1,1,1],
[1,0,1,0,1,1,1,0],
[0,1,0,1,0,1,1,1],
[1,1,1,1,0,0,1,0],
[0,1,1,1,1,0,0,1],
[1,1,1,0,0,1,0,1],
[1,0,1,0,1,0,1,1],
[1,0,0,0,1,1,0,0],
[0,1,0,0,0,1,1,0],
[0,0,1,0,0,0,1,1],
[1,1,0,0,1,0,0,0],
[0,1,1,0,0,1,0,0],
[0,0,1,1,0,0,1,0],
[0,0,0,1,1,0,0,1],
[1,1,0,1,0,1,0,1],
[1,0,1,1,0,0,1,1]])

tv1_tp4 = np.array([0,1,0,1,1,0,0,0,1,1,0,0,1,1,0,0,1,1,0,1,0,1,1,0,1,0,1,0,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,0,1,1,1,0,1,0,1,1,1,0,1,1,0,0,1,1,0,1,1,1,1,0,1,0,1,0,1,0,1,1,1,1,1,0,1,1,0,1,1,0,1,0,0,1,0,1,1,1,0,1,0,1,1,0,0,0,0,0,1,0,1,1,0,1,0,1,0,0,1,0,1,1,0,0,0,0]);
tv1_doubled = np.vstack((tv1_tp4, tv1_tp4))
tv1_parity = np.array([0,1,0,1,0,0,0,1])
eye_8 = np.eye(8, dtype = np.int32)
eye_60 = np.eye(60, dtype = np.int32)

#This is the generator matrix per https://www.ieee802.org/3/df/public/22_10/22_1005/bliss_3df_01_220929.pdf
G1 = np.vstack((eye_60, G.transpose()))
H1 = np.hstack((G.transpose(), eye_8))
def encode_177_5(G, M):
    # M is assumed to be a (2K) X 120
    
    ### From the draft:
    # The encode process is illustrated in Figure 177–5. Starting from the first message bit, each pair of bits that 
    # will be modulated into a PAM4 symbol are first XORed, resulting 60 bits, m_xor<59:0>, from the 120 bits 
    # of message. A Hamming encoder (68,60) will be used to generate the 8 parity bits based on the 60 bits 
    # obtained in the previous step.
    # For k = 0:59
    # m_xor(k) = m(2k) + m(2k+1)
    
    M_even = M[: , 0 :: 2]
    M_odd = M[: , 1 :: 2]
    M_xor = (M_even + M_odd) % 2
    M_parity = M_xor.dot(G) %2
    codewords = np.hstack(( M, M_parity))
    M_xor = np.hstack(( M_xor, M_parity))
    return codewords, M_xor

def encodeUsingMatrixOnly(G, M):
    return (G.dot(M) % 2)


def checkEqual():
    cw, mxor = encode_177_5(G, tv1_doubled)
    cw2 = encodeUsingMatrixOnly(G1, mxor[1,0:60].transpose())
    return np.all(cw2 == mxor[1,:].transpose())