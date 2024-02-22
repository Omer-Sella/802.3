# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 10:39:30 2024

@author: Omer
"""
import os
import sys
projectDir = os.environ.get('IEEE8032DJ')
if projectDir == None:
     projectDir = "D:/802.3/"
sys.path.insert(1, projectDir)
import numpy as np
from ieee8023dj_d0p1 import *


def simpleAWGNComparison():
    randomData = np.random.randint(0, 1, size = 1000)
    


