# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 15:30:10 2024

@author: Omer
"""

import numpy as np
import ieeeConstants
from ieee8023dj_d0p1 import encode_177_5
def test_encode_177_5_tv_1():
    tv1Encoded = encode_177_5(G = ieeeConstants.G, M = ieeeConstants.tv1_doubled)
    return tv1Encoded
    
    
