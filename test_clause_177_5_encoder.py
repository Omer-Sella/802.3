# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 15:30:10 2024

@author: Omer
"""

import numpy as np
import ieeeConstants
from ieee8023dj_d0p1 import encode_177_5
def test_encode_177_5_tv_1():
    tv1Encoded, tv1xored = encode_177_5(G = ieeeConstants.G, M = ieeeConstants.tv1_doubled)
    assert np.all(tv1Encoded[0,120:128] == ieeeConstants.tv1_parity)