# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 18:48:59 2023

@author: Megatron
"""
import numpy as np
import matplotlib.pyplot as plt
from channelFunctions import *
from modulationFunctions import *

SAMPLE_SIZE = 1000000
LOCAL_PRNG = np.random.RandomState(7134066)

def plotSNRvBER(SNRaxis, BERdata, fileName = None, inputLabel = 'baseline', figureNumber = 1, figureName = ''):
    snrBaseline = np.array([ 2. ,  2.5,  3. ,  3.5,  4. ,  4.5,  5. ,  5.5,  6. ,  6.5,  7. ,
         7.5,  8. ,  8.5,  9. ,  9.5, 10. ])
    berPam2 = np.array([3.75061284e-02, 2.96552876e-02, 2.28784076e-02, 1.71725417e-02,
        1.25008180e-02, 8.79381053e-03, 5.95386715e-03, 3.86223164e-03,
        2.38829078e-03, 1.39980484e-03, 7.72674815e-04, 3.98796335e-04,
        1.90907774e-04, 8.39995392e-05, 3.36272284e-05, 1.21088933e-05,
        3.87210822e-06])
    
    plt.ylabel('Output Bit Error Ratio (BER)',fontsize=16)
    plt.xlabel('Signal to noise ratio (SNR)',fontsize=16)
    plt.title(figureName)
    
    plt.semilogy(snrBaseline, berPam2, '--b', linewidth = 2, label = 'Uncoded PAM 2')
    #plt.plot(snrBaseline, berPam2, '--b', linewidth = 2, label = 'Uncoded PAM 2')
    
    plt.semilogy(SNRaxis, BERdata, '^', linewidth = 3, label=inputLabel)
    #plt.plot(SNRaxis, BERdata, '^', linewidth = 3, label=inputLabel)
    
    #plt.semilogy(snrActualNearEarthAxis, berNearEarthAxis, '*', linewidth = 3, label = 'NearEarthActual')
    #plt.plot(snrActualNearEarthAxis, berNearEarthAxis, '*', linewidth = 3, label = 'NearEarthActual')
    
    plt.tick_params(axis='both',  labelsize=16)
    #fig.set_size_inches(6.25, 6)
    plt.grid(True, which="both")
    plt.tight_layout()
    #plt.show()
    if fileName != None:
        plt.savefig(fileName, format='png', dpi=2000)
    #return timeStamp


def test_uncoded():
    SNRaxis = [ 2. ,  2.5,  3. ,  3.5,  4. ,  4.5,  5. ,  5.5,  6. ,  6.5,  7. ,
             7.5,  8. ,  8.5,  9. ,  9.5, 10. ]
    data = LOCAL_PRNG.randint(low = 0, high = 1, size = SAMPLE_SIZE)
    modulatedData = modulate(data, SAMPLE_SIZE)
    berData = np.array([])
    for s in SNRaxis:
        noisyData, _, _ = addAWGN(modulatedData, SAMPLE_SIZE, s, LOCAL_PRNG)
        slicedData = slicer(noisyData, SAMPLE_SIZE)
        ber = np.sum(slicedData != data) / SAMPLE_SIZE
        print(ber)
        berData = np.hstack((berData, ber))
    plotSNRvBER(SNRaxis, berData)
    return 'OK'