# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 10:56:49 2019

@author: mugdhapolimera

Plotting spectrum from SDSS MPA-JHU table values
"""
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import norm
import numpy as np
from scipy.interpolate import interp1d
import os

os.chdir('C:\Users\mugdhapolimera\github\SDSS_spectra')

df = pd.read_csv('RESOLVE_full_cont.csv',index_col = 'name')
galname = 'rs0010'
gal = df.loc[galname]

wavelengths = {'oii_3726' : 3726.032, 'oii_3729' : 3728.815, 
               'neiii_3869' : 3868.760, 'h_gamma' : 4340.471, 
               'oiii_4363' : 4363.210, 'heii_4685' : 4685.710,
               'h_beta' : 4861.333, 'oiii_4959': 4958.911, 
               'oiii_5007' : 5006.843, 'hei_5875' : 5875.624, 
               'oi_6300' : 6300.304, 'nii_6548' : 6548.050, 
               'h_alpha' : 6562.819, 'nii_6584' : 6583.46, 
               'sii_6717' : 6716.440, 'sii_6731' : 6730.810,
               'ariii_7136' : 7135.790}

newlam = np.arange(3724,7400,0.5)
emspec = np.zeros(len(newlam))
cont = [] #np.zeros(len(newlam))
#plt.figure()
for lam in wavelengths.keys():
    mu = wavelengths[lam]
    if lam+'_flux' in gal.keys():
        sigma = 1.69/2.355 #FWHM to sigma (69/3e5)*6500
        x = (newlam - mu)/sigma
        #including flux conservation as SDSS_FWHM/GEMINI_FWHM,
        line = norm.pdf(x)#*np.array(gal[lam+'_flux'])#*(3.6/1.69) 
        #print line, max(line), lam
        if np.max(line):
            line = (line/np.max(line))*0.4*np.array(gal[lam+'_flux'])/(sigma)#*(3.6/1.69)         
        emspec+= line
        cont.append([mu,gal[lam+'_cont']])
        #plt.plot(newlam,line,'b')
cont = np.array(cont)
cont = cont[cont[:,0].argsort()]
contnew = interp1d(cont[:,0],cont[:,1], fill_value = 'extrapolate')(newlam)
#plt.plot(cont[:,0],cont[:,1],'r')
#plt.plot(newlam,contnew,'g')

emspec += contnew
plt.figure()
plt.plot(newlam,emspec*1e-17,'r')
plt.ylabel('Flux ergs/s/cm^2')
plt.xlabel('Wavelength (Angstroms)')
optndx = (newlam < 7400) & (newlam > 3000)
os.chdir('C:\Users\mugdhapolimera\github\BPT\grids\gemini\/')
np.savetxt('mock_SDSStable_rs0010.txt', \
           zip(newlam[optndx]/10,emspec[optndx]*10*1e-17))
