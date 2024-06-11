# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 09:44:44 2024

@author: Omer
"""
import numpy as np
from ieeeConstants import *
def encode_177_5(G, M):
    
    # The encode process is illustrated in Figure 177–5. Starting from the first message bit, each pair of bits that 
    # will be modulated into a PAM4 symbol are first XORed, resulting 60 bits, m_xor<59:0>, from the 120 bits 
    # of message. A Hamming encoder (68,60) will be used to generate the 8 parity bits based on the 60 bits 
    # obtained in the previous step.
    # For k = 0:59
    # m_xor(k) = m(2k) + m(2k+1)
    
    M_xor = np.array([ (M[2*k] + M[2 * k + 1]) %2 for k in range(60) ]) # Hardcoded - the length of the data message is 120 !
    M_parity = M_xor.dot(G) %2
    M_xor = np.expand_dims(M_xor, axis = 0)
    print(M_parity.shape)
    print(M_xor.shape)
    codeword = np.hstack(( M_xor, M_parity))
    return codeword




def bchEncoder(M, G = None):
    #M is assumed to be a matrix of size L X 110 if G is not None, and where if G is None it is assumed that L == 1
    #If G is not none then it is assumed to be a binary generator matrix
    """
     184.4.5 BCH encoder
     The BCH encoder shall work in conjunction with the outer RS(544,514) FEC to provide a high-performance 
    FEC for 800GBASE-LR1. There are 32 BCH encoder functions.
     The BCH code defined by the generator polynomial is derived starting with a BCH(127,113) code that is 
    shortened to a BCH(124,110) code and then extended to the final BCH(126,110) code by adding even and 
    odd parity check bits.
     The code is defined as follows:
     — Define c as a binary vector of length 126, and c(x) a polynomial of degree 125 with the coefficients 
    defined by c (where the bit 0 of c represents the coefficient of power 125).
     — Then c is a codeword of the BCH(126,110) code if c(x) is divisible (modulo 2) by the defined binary 
    generator polynomial g(x)
     BCH (126,110) code operates on 11 RS symbols provided by the convolutional interleaver
     — Denote the index of the aligned 40-bit block at the convolutional interleaver output as i and the bit 
    index within block i as j where j=0 to 39. The indices i and j are identical to the indices used in the 
    definition of the convolutional interleaver.
     — Denote the index of the BCH encoder block as u and the bit index within the input block as v where 
    v=0 to 109. The indices u and v are related to i and j and defined as u = [(40i + j) / 110] and 
    v = mod((40i + j), 110).
     The encoding of each block u on lane p is defined as follows:
     — A message polynomial m(x) of degree 109 is defined where the coefficient of the x109-v term is the 
    BCH encoder input_p [110 u + v] for v = 0 to 109 
    — A generator polynomial g(x) of degree 16 is defined as 
    g(x) = x16 +x 14 + x11 + x10 + x9 + x7 + x5 + x3 + x+ 1
     — A parity polynomial p(x) of degree 15 is defined as the remainder from the division (module 2) of 
    m(x) x16 by the generator polynomial g(x) . p(x) =p15x15 + p14x14 + .... +p1x + .p0
     The bit mapping of the BCH encoder output is:
     For v=0 to 109
     output_p[126u + v] = input_p[110u + v] 
    For v=110 to 125
     output_p[126u + v] = p125-v
    """
    if G is not None:
        codeword = G.dot(M) % 2
    else:
        # g(x) =      x16 +  x14 +    x11 + x10 + x9 + x7 + x5 +  x3 +  x+ 1
        gX = np.array([1, 0, 1, 0, 0, 1,    1,    1,0 ,1, 0,1, 0, 1, 0, 1, 1]) 
        remainder = np.zeros((len(M) + len(gX) - 1), dtype = IEEE_8023_INT_DATA_TYPE) # Padded with 15 zeros
        #print(len(M))
        remainder[0:len(M)] = M
        #print(remainder)
        for i in range(len(M)):
            #print(i)
            if remainder[i] == 1:
                #kill ther leading coefficient
                #print('killing the leader at coordinate:' + str(i))
                #print((remainder[i : i + (len(gX))] + gX) % 2)
                remainder[i : i + len(gX)] = (remainder[i : i + (len(gX))] + gX) % 2
        codeword = np.zeros((len(M) + len(gX) - 1), dtype = IEEE_8023_INT_DATA_TYPE)
        codeword[0 : len(M)] = M
        codeword[len(M) : len(codeword)] = remainder[len(remainder) - len(gX) + 1 : ]
    return codeword
        


