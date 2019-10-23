# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 09:29:21 2019

@author: mugdhapolimera

Plot model spectra from Cloudy for Z = 0.4 Z_solar, log q = 7.9, M_bh = 10^5 M_sun
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from PyAstronomy import pyasl
from PyAstronomy import funcFit as fuf

from pysynphot import observation
from pysynphot import spectrum
from scipy.io import readsav
 
def rebin_spec(wave, specin, wavnew):
    spec = spectrum.ArraySourceSpectrum(wave=wave, flux=specin)
    f = np.ones(len(wave))
    filt = spectrum.ArraySpectralElement(wave, f, waveunits='angstrom')
    obs = observation.Observation(spec, filt, binset=wavnew)
 
    return obs.binflux

#agnfrac = [0, 4, 16, 32, 50, 64, 100]
#spec = pd.read_csv('agn_0000-Z_0400.con', sep = '\t')
#spec = spec.rename(columns = {'#Cont  nu':'lam'})
#plt.plot(spec.lam, spec.total)
#plt.plot(6545.97*np.ones(len(np.linspace(0,max(spec.total),10))),(np.linspace(0,max(spec.total),10)),'r-.')
#plt.plot(6589.76*np.ones(len(np.linspace(0,max(spec.total),10))),(np.linspace(0,max(spec.total),10)),'r-.')
#plt.xlim(4000,8000)
#spectra-bin.z008.dat spectra-bin-imf135all_100.z008
#spec = pd.read_csv('spectra-bin.z008_1.dat',sep = '\s+', header=None, 
#                   usecols = [0,41], names = ['lam', 'flux'])
spec = pd.read_csv('agn_0000-Z_0400.con',sep = '\s+', header=None, skiprows = 1,
                   usecols = [0,6], names = ['lam', 'flux'])
spec = spec.iloc[np.where((spec.lam < 10000 ) & (spec.lam >=1 ))]
#spec = spec.iloc[np.where((8000 > spec.lam) & (spec.lam > 4000))]
plt.figure()
plt.plot(spec.lam,spec.flux/max(spec.flux))
newlam = np.arange(min(spec.lam),max(spec.lam),0.45)
spec.flux = np.array(spec.flux[::-1])
spec.lam = np.array(spec.lam[::-1])
new_spec = rebin_spec(np.array(spec.lam),np.array(spec.flux),newlam)
plt.plot(newlam,new_spec/max(new_spec), 'r')

#vacuum to air conversion
s = 10**4/newlam
n = 1 + 0.0000834254 + 0.02406147 / (130 - s**2) + 0.00015998 / (38.9 - s**2)
#airlam = newlam/n 
#air_spec = rebin_spec(np.array(newlam),np.array(new_spec),airlam)
#new_spec = air_spec
#plt.plot(airlam,new_spec/max(new_spec), 'k')
#new_spec, vind = pyasl.specAirVacConvert(newlam, new_spec, \
#                direction="vactoair")
plt.plot(newlam,new_spec/max(new_spec), 'k')

resdata = readsav('../../../SDSS_spectra/resolvecatalog.dat')
galname = 'rs0010'

v = resdata['vhel'][np.where(resdata['name'] == galname)] #km/s
z = v/3e5 #redshift
r_pc = (v/70)*10**6*(1+z) #pc
r = r_pc*3.086e+18 #pc to cm
mstar = resdata['mstars'][np.where(resdata['name'] == galname)]
#new_spec = new_spec/4*np.pi*r
emlines = pd.read_csv('../../../izi/Richardson-0-0_1-0agn-M-5_0-BPASS-Binary-CSF-n=1e3-40.0Myr-NichollsCE-D_G-RR14_Fstar_0_3-unified-2.csv')
emlines = emlines[(emlines.AGNFRAC == 0.0) & (emlines.LOGZ == np.log10(0.4)) &\
                   (7.3>emlines.LOGQ) & (emlines.LOGQ > 7.2)]
wavelengths = {'oii3726' : 3726.032, 'oii3729' : 3728.815, 
               'neiii3869' : 3868.760, 'hgamma' : 4340.471, 
               'oiii4363' : 4363.210, 'heii4685' : 4685.710,
               'hbeta' : 4861.333, 'oiii4959': 4958.911, 
               'oiii5007' : 5006.843, 'hei5875' : 5875.624, 
               'oi6300' : 6300.304, 'nii6548' : 6548.050, 
               'halpha' : 6562.819, 'nii6584' : 6583.46, 
               'sii6717' : 6716.440, 'sii6731' : 6730.810,
               'ariii7136' : 7135.790}

emspec = np.zeros(len(new_spec))
#plt.figure()
for lam in wavelengths.keys():
    mu = wavelengths[lam]
    if mu > newlam[0]:
        sigma = 1.69/2.355 #FWHM to sigma
        x = (newlam - mu)/sigma
        line = norm.pdf(x)*np.array(emlines[lam])
        #print line, max(line), lam
        line = (line/max(line)) * np.array(emlines[lam]) 
        emspec+= line
        plt.plot(newlam,line,'g')

#totalspec = (new_spec/max(new_spec))+(emspec/max(emspec))
totalspec = new_spec*(1+emspec)
plt.plot(newlam,totalspec,'m')
lsun = 3.826*10**33 #ergs/s
age = 10**(6+0.1*(42-2)) # in years
norm =1# 25*1e-2
#totalflux = norm*mstar*(totalspec* lsun) /(4*np.pi*r*r*1e6)
totalflux = new_spec/(newlam)

plt.figure()
plt.plot(newlam, totalflux)
plt.xlabel('Wavelength (in Angstroms)')
plt.ylabel('Specific Flux (in ergs/s/cm^2/Angstrom)')
#plt.xlim(0,10000)
#add reddening
#Assume E_BV of one of our galaxies?

#add broadening

