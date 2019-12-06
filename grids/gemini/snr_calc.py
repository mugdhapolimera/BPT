# -*- coding: utf-8 -*-
"""
Created on Sat Oct 26 19:29:13 2019

@author: mugdhapolimera

Calculate the S/N required to be able to measure line broadening to a certain
confidence, given the width of a Gaussian
"""
import numpy as np

#FWHM of the Broad component of Halpha ~ 500 km/s for LMBH
#Width (sigma) of the Broad Halpha line in Angstroms (observed)
#sigma = (velocity_FWHM/2.354)/c x wavelength
vel_FWHM = 500
sig_b = ((vel_FWHM/2.354)/3e5)*6563
#Variance,
sig_b2 = sig_b**2
#Width (sigma) of the Broad Halpha line in Angstroms (observed) ~ instrumental 
#resolution - GMOS IFU - 1.69 FWHM 
sig_n = 1.69/2.354
#We need to detect sig_b at a 5-sigma confidence 
#(sig_b - sig_n)/sigma(sig_b) >= 5
confidence = 5
sig_sig_b = (sig_b - sig_n)/confidence
sig_sig_b = sig_b/confidence
sig_sig_b = ((33)/3e5)*6563
#This means that sig_b +- sig_sig_b is needed to get a 5-sigma confirmation
#of a 500 km/s FWHm of broad Halpha component

#sig_sig_b is basically error(sig_b). To convert this to variance, we use 
#error propagation to get error(sig_b^2) = 2 * sig_b * error(sig_b)
sig_sig_b2 = 2*sig_b*sig_sig_b

#Variance of the variance of broad Halpha line,
var_sig_b2 = sig_sig_b2**2

#Sample variance, sb^2 = sig_b**2 * (N-1)/N
#Variance of the sample variance, var(s_b2) = (N-1/N)**2 var(sig_b2) 
#Variance of the sample variance is related to sigma of the distribution by,
#var(s_b2) =  2 sig_b**4 * (N-1)/N**2

#(N-1/N)**2 var(sig_b2) = 2 sig_b**4 * (N-1)/N**2
#(N-1) var(sig_b2) = 2 sig_b**4
N = 1 + 2*sig_b**4/var_sig_b2
SNR = np.ceil(np.sqrt(N))
print(np.sqrt(N))
print('For broad Halpha component with FWHM velocity {} km/s, we need S/N >= {}'.format(vel_FWHM, SNR, confidence))# for a {}-sigma detection'
print('1-sigma error on the width of broad Halpha is {}km/s'.format(sig_sig_b*3e5/6563))