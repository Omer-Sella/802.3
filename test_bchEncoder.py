# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 14:26:43 2024

@author: Omer
"""

import numpy as np
from ieeeConstants import *
from ieee8023dj_d0p1 import bchEncoder
def test_bchEncoder_zero_data():
    zeroData = np.zeros(110, dtype = IEEE_8023_INT_DATA_TYPE)
    #tv0[109] = 1
    zeroDataEncoded = bchEncoder(zeroData)
    assert (np.all(zeroDataEncoded == 0))
    
    
