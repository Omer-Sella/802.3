

import scipy.io
import sys
import os
import numpy as np

projectDir = os.environ.get('8023')
if projectDir == None:
     projectDir = "D:/802.3/"
sys.path.insert(1, projectDir)

workspace = scipy.io.loadmat('./bliss_3df_01_220929.mat')

P = workspace['p']

I_120_120 = np.identity(120)
G = np.vstack((P,I_120_120))

H = np.hstack((np.identity(8), P))

# Check H is a parity matrix for G:
assert(np.all(H.dot(G) % 2 ==0))