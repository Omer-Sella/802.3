# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 14:37:58 2024

@author: Omer
"""
import matplotlib.pyplot as plt
import numpy as np

def plotSNRvsBER(SNRaxis, BERdata, fileName = None, inputLabel = 'baselineCode', figureNumber = 1, figureName = '', figureHandle = None):
    snrBaseline = np.array([ 2. ,  2.5,  3. ,  3.5,  4. ,  4.5,  5. ,  5.5,  6. ,  6.5,  7. ,
         7.5,  8. ,  8.5,  9. ,  9.5, 10. ])
    berPam2 = np.array([3.75061284e-02, 2.96552876e-02, 2.28784076e-02, 1.71725417e-02,
        1.25008180e-02, 8.79381053e-03, 5.95386715e-03, 3.86223164e-03,
        2.38829078e-03, 1.39980484e-03, 7.72674815e-04, 3.98796335e-04,
        1.90907774e-04, 8.39995392e-05, 3.36272284e-05, 1.21088933e-05,
        3.87210822e-06])
    #timeStamp = str(time.time())
    #fileName = "./" + timeStamp + fileName
    if figureHandle == None:
        fig, ax = plt.subplots()
        #plt.semilogy(snrBaseline, berPam2, '--b', linewidth = 2, label = 'Uncoded PAM 2')
        fig.plot(snrBaseline, berPam2, '--b', linewidth = 2, label = 'Uncoded PAM 2')
        fig.ylabel('Output Bit Error Ratio (BER)',fontsize=16)
        fig.xlabel('Signal to noise ratio (SNR)',fontsize=16)
        fig.title(figureName)    
        fig.tick_params(axis='both',  labelsize=16)
        #fig.set_size_inches(6.25, 6)
        fig.grid(True, which="both")
        #fig.tight_layout()
    
    
    #plt.semilogy(SNRaxis, BERdata, '^', linewidth = 3, label=inputLabel)
    fig.plot(SNRaxis, BERdata, '^', linewidth = 3, label=inputLabel)
       
    #plt.semilogy(snrActualNearEarthAxis, berNearEarthAxis, '*', linewidth = 3, label = 'NearEarthActual')
    #plt.plot(snrActualNearEarthAxis, berNearEarthAxis, '*', linewidth = 3, label = 'NearEarthActual')
    
    
    #plt.show()
    if fileName != None:
        plt.savefig(fileName, format='png', dpi=2000)
    return fig, ax
    