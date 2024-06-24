# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 20:26:57 2024

@author: Omer
"""

import scipy.io
import sys
import os
projectDir = os.environ.get('IEEE8023DJ')
if projectDir == None:
     projectDir = "c:/users/omer/802.3/"
sys.path.insert(0, projectDir)

import numpy as np
IEEE_8023_DATA_TYPE = np.int64
IEEE_8023_INT_DATA_TYPE = np.int64
IEEE_8023_DECIMAL_DATA_TYPE = np.float64
IEEE_8023_DATA_TYPE = np.int64

# Omer Sella: seeds can be integers between 0 and 2**31 - 1
IEEE_8023_MAX_SEED = 2**31 - 1

#NUMBA_INT = int64
#NUMBA_FLOAT = float64
#NUMBA_BOOL = boolean

PAM4_LEVEL_LOW = -1
PAM4_LEVEL_MID_LOW = -(1/3)
PAM4_LEVEL_MID_HIGH = 1/3
PAM4_LEVEL_HIGH = 1

PAM4_LEVELS = np.array([PAM4_LEVEL_LOW, PAM4_LEVEL_MID_LOW, PAM4_LEVEL_MID_HIGH, PAM4_LEVEL_HIGH])

PAM4_GRAYCODED = np.array([[0,0],[0,1],[1,1],[1,0]])
PAM4_NOGRAYCODING = np.array([[0,0],[0,1],[1,0],[1,1]])
                          
              

g_177_1 = np.array(
  #Table 177–1—Generation matrix for Hamming(68,60)
 [
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
tv1_parity = np.array([0,1,0,1,0,0,0,1])
eye_8 = np.eye(8, dtype = np.int32)
eye_60 = np.eye(60, dtype = np.int32)

#This is the generator matrix per https://www.ieee802.org/3/df/public/22_10/22_1005/bliss_3df_01_220929.pdf
# Make a 68,60 Hamming generator matrix from the 60,8 matrix in the draft
generatorMatrix_177_5 = np.vstack((eye_60, g_177_1.transpose()))
# Now create a parity matrix using the instructions from https://www.ieee802.org/3/df/public/22_10/22_1005/bliss_3df_01_220929.pdf
parityMatrix_177_5 = np.hstack((g_177_1.transpose(), eye_8))
assert(np.all(parityMatrix_177_5.dot(generatorMatrix_177_5) % 2 ==0))

snrBaseline = np.array([ 1. ,  1.2,  1.4,  1.6,  1.8,  2. ,  2.2,  2.4,  2.6,  2.8,  3. ,
             3.2,  3.4,  3.6,  3.8,  4. ,  4.2,  4.4,  4.6,  4.8,  5. ,  5.2,
             5.4,  5.6,  5.8,  6. ,  6.2,  6.4,  6.6,  6.8,  7. ,  7.2,  7.4,
             7.6,  7.8,  8. ,  8.2,  8.4,  8.6,  8.8,  9. ,  9.2,  9.4,  9.6,
             9.8, 10. , 10.2, 10.4, 10.6, 10.8, 11. , 11.2, 11.4, 11.6, 11.8,
            12. , 12.2, 12.4, 12.6, 12.8, 13. , 13.2, 13.4, 13.6, 13.8, 14. ,
            14.2, 14.4, 14.6, 14.8, 15. , 15.2, 15.4, 15.6, 15.8, 16. ])

berPam2 = np.array([5.6282e-02, 5.2216e-02, 4.8301e-02, 4.4541e-02, 4.0942e-02,
            3.7506e-02, 3.4238e-02, 3.1140e-02, 2.8214e-02, 2.5460e-02,
            2.2878e-02, 2.0469e-02, 1.8229e-02, 1.6157e-02, 1.4249e-02,
            1.2501e-02, 1.0907e-02, 9.4624e-03, 8.1600e-03, 6.9930e-03,
            5.9539e-03, 5.0346e-03, 4.2269e-03, 3.5223e-03, 2.9123e-03,
            2.3883e-03, 1.9419e-03, 1.5648e-03, 1.2492e-03, 9.8751e-04,
            7.7267e-04, 5.9812e-04, 4.5782e-04, 3.4634e-04, 2.5880e-04,
            1.9091e-04, 1.3894e-04, 9.9706e-05, 7.0501e-05, 4.9086e-05,
            3.3627e-05, 2.2650e-05, 1.4989e-05, 9.7362e-06, 6.2027e-06,
            3.8721e-06, 2.3663e-06, 1.4142e-06, 8.2572e-07, 4.7048e-07,
            2.6131e-07, 1.4130e-07, 7.4294e-08, 3.7934e-08, 1.8783e-08,
            9.0060e-09, 4.1752e-09, 1.8686e-09, 8.0599e-10, 3.3447e-10,
            1.3329e-10, 5.0917e-11, 1.8606e-11, 6.4904e-12, 2.1566e-12,
            6.8102e-13, 2.0389e-13, 5.7727e-14, 1.5416e-14, 3.8725e-15,
            9.1240e-16, 2.0101e-16, 4.1281e-17, 7.8764e-18, 1.3914e-18,
            2.2674e-19])
berPam2Hamming_127_120 = np.array([6.1465e-02, 5.7208e-02, 5.3088e-02, 4.9107e-02, 4.5266e-02,
        4.1566e-02, 3.8006e-02, 3.4584e-02, 3.1301e-02, 2.8155e-02,
        2.5147e-02, 2.2281e-02, 1.9564e-02, 1.7003e-02, 1.4610e-02,
        1.2398e-02, 1.0377e-02, 8.5590e-03, 6.9488e-03, 5.5480e-03,
        4.3526e-03, 3.3529e-03, 2.5344e-03, 1.8787e-03, 1.3650e-03,
        9.7167e-04, 6.7735e-04, 4.6221e-04, 3.0860e-04, 2.0152e-04,
        1.2865e-04, 8.0241e-05, 4.8873e-05, 2.9050e-05, 1.6839e-05,
        9.5116e-06, 5.2307e-06, 2.7979e-06, 1.4542e-06, 7.3356e-07,
        3.5871e-07, 1.6981e-07, 7.7717e-08, 3.4335e-08, 1.4620e-08,
        5.9896e-09, 2.3570e-09, 8.8918e-10, 3.2097e-10, 1.1063e-10,
        3.6333e-11, 1.1343e-11, 3.3588e-12, 9.4090e-13, 2.4871e-13,
        6.1864e-14, 1.4439e-14, 3.1530e-15, 6.4210e-16, 1.2155e-16,
        2.1316e-17, 3.4505e-18, 5.1363e-19, 7.0031e-20, 8.7101e-21,
        9.8392e-22, 1.0049e-22, 9.2375e-24, 7.6055e-25, 5.5599e-26,
        3.5988e-27, 2.1516e-28, 1.2030e-29, 8.3816e-31, 0.0000e+00,
        0.0000e+00])
berPam4 = np.array([1.1900e-01, 1.1468e-01, 1.1040e-01, 1.0615e-01, 1.0193e-01,
            9.7742e-02, 9.3594e-02, 8.9488e-02, 8.5426e-02, 8.1413e-02,
            7.7453e-02, 7.3551e-02, 6.9712e-02, 6.5941e-02, 6.2243e-02,
            5.8624e-02, 5.5089e-02, 5.1643e-02, 4.8291e-02, 4.5040e-02,
            4.1893e-02, 3.8855e-02, 3.5930e-02, 3.3122e-02, 3.0435e-02,
            2.7871e-02, 2.5433e-02, 2.3123e-02, 2.0941e-02, 1.8889e-02,
            1.6967e-02, 1.5173e-02, 1.3506e-02, 1.1965e-02, 1.0546e-02,
            9.2472e-03, 8.0637e-03, 6.9913e-03, 6.0253e-03, 5.1602e-03,
            4.3903e-03, 3.7098e-03, 3.1123e-03, 2.5914e-03, 2.1409e-03,
            1.7542e-03, 1.4250e-03, 1.1472e-03, 9.1491e-04, 7.2250e-04,
            5.6471e-04, 4.3664e-04, 3.3382e-04, 2.5222e-04, 1.8822e-04,
            1.3866e-04, 1.0077e-04, 7.2207e-05, 5.0977e-05, 3.5435e-05,
            2.4234e-05, 1.6294e-05, 1.0762e-05, 6.9771e-06, 4.4359e-06,
            2.7632e-06, 1.6848e-06, 1.0046e-06, 5.8510e-07, 3.3252e-07,
            1.8419e-07, 9.9315e-08, 5.2065e-08, 2.6502e-08, 1.3080e-08,
            6.2502e-09])
berPam4Hamming_127_120 = np.array([1.2435e-01, 1.2000e-01, 1.1568e-01, 1.1139e-01, 1.0713e-01,
        1.0290e-01, 9.8706e-02, 9.4549e-02, 9.0432e-02, 8.6360e-02,
        8.2334e-02, 7.8361e-02, 7.4443e-02, 7.0586e-02, 6.6795e-02,
        6.3074e-02, 5.9427e-02, 5.5858e-02, 5.2372e-02, 4.8971e-02,
        4.5657e-02, 4.2434e-02, 3.9302e-02, 3.6262e-02, 3.3315e-02,
        3.0462e-02, 2.7704e-02, 2.5044e-02, 2.2486e-02, 2.0036e-02,
        1.7703e-02, 1.5494e-02, 1.3421e-02, 1.1495e-02, 9.7251e-03,
        8.1203e-03, 6.6858e-03, 5.4237e-03, 4.3318e-03, 3.4038e-03,
        2.6297e-03, 1.9964e-03, 1.4885e-03, 1.0894e-03, 7.8232e-04,
        5.5095e-04, 3.8035e-04, 2.5727e-04, 1.7042e-04, 1.1050e-04,
        7.0090e-05, 4.3467e-05, 2.6337e-05, 1.5580e-05, 8.9917e-06,
        5.0581e-06, 2.7708e-06, 1.4766e-06, 7.6463e-07, 3.8432e-07,
        1.8725e-07, 8.8319e-08, 4.0270e-08, 1.7722e-08, 7.5162e-09,
        3.0667e-09, 1.2016e-09, 4.5128e-10, 1.6214e-10, 5.5614e-11,
        1.8171e-11, 5.6430e-12, 1.6616e-12, 4.6274e-13, 1.2157e-13,
        3.0045e-14])

#  120.5.11.2.a PRBS9Q test pattern
# For example, if the PRBS9 generator used to create the PRBS9Q
# sequence is initialized to a seed value of 111111111 (with the leftmost bit in S0 and the rightmost in S8), the 
# PRBS9Q sequence is the following sequence of Gray coded PAM4 symbols, transmitted left to right:
stringValuePRBS9Q_seed111111111_clause_120_5_11_2_a = "0012322303231310010331213302202231320111030230213332303130303000100302003120333200212331323101100332102221310311322203133313130002013110133112221011302332032022012212100133213232001133223333300110332203232300120233102211211010301312003221320210023220022223"
stringValuePRBS9Q_seed111111111_clause_120_5_11_2_a += "002212201120203003110232101231220213033310120132111201020101000030101301023111130132210212030330111331223203103212231021102020001302033021032223303201211311312302232330021132121300321122111100033111231121200023121031233233303100202301123213133012123012222"
PRBS9Q_seed111111111_clause_120_5_11_2_a = [int(i) for i in list(stringValuePRBS9Q_seed111111111_clause_120_5_11_2_a)]
