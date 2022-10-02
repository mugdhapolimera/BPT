# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 18:38:18 2022

@author: mugdhapolimera
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 09:29:21 2019

@author: mugdhapolimera

Plot model spectra from Cloudy for Z = 0.4 Z_solar, log q = 7.9, M_bh = 10^5 M_sun
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
from scipy.stats import norm
from PyAstronomy import pyasl
from PyAstronomy import funcFit as fuf
from scipy import integrate

from pysynphot import observation
from pysynphot import spectrum
from scipy.io import readsav

def stern(x):
    return 0.8*np.ones(len(x))

def jarretx(y):
    return [2.2*np.ones(len(y)), 4.2*np.ones(len(y))]

def jarrety(x):
    return [1.7*np.ones(len(x)), 0.1*x+0.38]

def satyapalx(x):
    return 0.52*np.ones(len(x))

def satyapaly(x):
    return 5.78*x -24.50

#re-bin the spectra while conserving flux 
def rebin_spec(wave, specin, wavnew):
    spec = spectrum.ArraySourceSpectrum(wave=wave, flux=specin)
    f = np.ones(len(wave))
    filt = spectrum.ArraySpectralElement(wave, f, waveunits='angstrom')
    obs = observation.Observation(spec, filt, binset=wavnew)
 
    return obs.binflux

spectype = 'cloudy'

folder = r'C:\Users\mugdhapolimera\github\SDSS_spectra\mid_ir\satyapalmodels\/'

#Specific models 
#speclowmbh = pd.read_csv(folder+r'grid000000002_Z_1_n_300_100per.con',sep = '\s+', header=None, skiprows = -2,
#                   usecols = [0,6], names = ['lam', 'flux'])[:-1] #agn_0000-Z_0400.con
#speclowmbh = speclowmbh.astype(np.float64)
#spechighmbh = pd.read_csv(folder+r'grid000000013_Z_1_n_300_100per.con',sep = '\s+', header=None, skiprows = 1,
#               usecols = [0,6], names = ['lam', 'flux'])[:-1] #agn_0000-Z_0400.con
#spechighmbh = spechighmbh.astype(np.float64)

#All models
bhspecfiles = glob.glob(folder+'*.con')

#Plot spectra as lam*F_lam vs. lam (in microns)
figspec, specax = plt.subplots(1)
specax.set_xlabel(r'Wavelength (in ${\mu m}$)')
specax.set_ylabel(r'Flux Density ')#($10^{-17}$ $ergs/s/cm^2/{\AA}$)')
specax.set_yscale('log')
specax.set_xscale('log')
specax.set_xlim(0.1,100)
#    specax.ylim(0.001,10)
specax.axvspan(3,3.8, alpha = 0.3, color = 'k')
specax.axvspan(4,5, alpha = 0.3, color = 'k')
specax.axvspan(9,13, alpha = 0.3, color = 'k')
specax.axvspan(20,25, alpha = 0.3, color = 'k')

#Plot mid ir colour cuts
figcol,ax = plt.subplots(1)
xaxis = np.linspace(0,6)
yaxis = np.linspace(jarrety(np.array([2.2]))[1],1.7)
#ax.plot(xaxis, stern(xaxis), 'k-.')#, label = 'Stern12')
#ax.text(5.75,1.0,'St12', fontsize = 15)
ax.plot(xaxis, satyapalx(xaxis), 'k')#, label = 'Satyapal18')
ax.text(5.75,0.3,'Sa14', fontsize = 15)
ax.text(4.7,2.0 ,'Sa18', fontsize = 15)
#xaxis = np.linspace(4.3287,6)
#ax.plot(xaxis, satyapaly(xaxis), 'k')

xaxis = np.linspace(2.2,4.2)
ax.plot(jarretx(yaxis)[0], yaxis, 'k--', jarretx(yaxis)[1], yaxis, 'k--')
ax.plot(xaxis, jarrety(xaxis)[0], 'k--')
ax.plot(xaxis, jarrety(xaxis)[1],'k--')#, label = 'Jarrett15')
ax.text(3.5,1.85,'J11', fontsize = 15)
ax.set_xlabel('W2 - W3')
ax.set_ylabel('W1 - W2')
ax.set_ylim(0, 2.5)
ax.set_xlim(1,6.5)

for specfile in bhspecfiles:
    spec = pd.read_csv(specfile,sep = '\s+', header=None, skiprows = -2,
                       usecols = [0,6], names = ['lam', 'flux'])[:-1] 
    spec = spec.astype(np.float64)
    spec = spec.iloc[np.where((spec.lam < 1e6 ) & (spec.lam >=1000 ))]

    #make finer spacing for the spectra
    newlam = np.arange(np.min(spec.lam),np.max(spec.lam),0.5)
    spec.flux = np.array(spec.flux[::-1])
    spec.lam = np.array(spec.lam[::-1])
    new_spec = rebin_spec(np.array(spec.lam),np.array(spec.flux),newlam)
    
    #vacuum to air conversion
    s = 10**4/newlam
    n = 1 + 0.0000834254 + 0.02406147 / (130 - s**2) + 0.00015998 / (38.9 - s**2)
    airlam = newlam/n 
    new_spec = new_spec*n
    newlam = airlam
    
#   If needed, normalize by hbeta. ignoring this by default    
#    hbeta_flux = new_spec[np.where(newlam == newlam[newlam>=4864.55][0])]
#    new_spec = new_spec/hbeta_flux
#    plt.plot(newlam,new_spec, 'k')

#Read SF-AGN recession velocity from the database - ignoring for now and hardcoding 
#the velocity below
#    resdata = readsav('../../../SDSS_spectra/resolvecatalog.dat')
#    galname = 'rs0010'
#    resphot = readsav('../../../SDSS_spectra/resolvecatalogphot.dat')
#    galndx = np.where(resdata['name'] == galname)
#    v = resdata['vhel'][galndx] #km/s

#hardcoded recession velocity of an example SF-AGN
    v = 5828.03 #km/s
    z = v/3e5 #redshift
    r_pc = (v/70)*1e6 *(1+z) #pc
    r = np.float64(r_pc*3.086e+18) #pc to cm

    #Convert lam*L_lam to lam*F_lam
    totalflux = new_spec/4*np.pi*r*r #units of ergs/s/cm^2 -- nu*F_nu or lam*F_lam
    totalflux = totalflux*1e-7 #convert to W/cm^2

    newlam = newlam*1e-10/1e-6 #convert angstrom to microns
    totalflux = totalflux/newlam #CONVERT LAM*f_LAM to f_lam

    #to redshift the spectrum - TBD; need to use GR to do this properly. ignore for now
#    totalflux = totalflux/(1+z) #lam = lam*(1+z); flux = flux/(1+z); flux density lam*F_lam stays same
#    newlam = newlam*(1+z)
    specax.plot(newlam, totalflux)

#defining the wavelength limits of each of the bands    
    w1ndx = (newlam > 3) & (newlam < 3.8)
    w2ndx = (newlam > 4) & (newlam < 5)
    w3ndx = (newlam > 9) & (newlam < 13)
    w4ndx = (newlam > 20) & (newlam < 25)

#read in response functions from file and convolve them with spectra
    w1response = pd.read_csv(r'w1response.txt', sep = '\s+')
    w1ndx = (newlam > np.min(w1response.W1)) & (newlam < np.max(w1response.W1))
    lamF_lam_w1 = rebin_spec(np.array(newlam[w1ndx]), \
                          np.array(totalflux[w1ndx]),w1response.W1)
    fluxmw1 = integrate.simps((w1response.RSR*lamF_lam_w1), w1response.W1)#/6.6e-1

    w2response = pd.read_csv(r'w2response.txt', sep = '\s+')
    w2ndx = (newlam > np.min(w2response.W2)) & (newlam < np.max(w2response.W2))
    lamF_lam_w2 = rebin_spec(np.array(newlam[w2ndx]), \
                          np.array(totalflux[w2ndx]),w2response.W2)
    fluxmw2 = integrate.simps((w2response.RSR*lamF_lam_w2), w2response.W2)#/1.0423

    w3response = pd.read_csv(r'w3response.txt', sep = '\s+')
    w3ndx = (newlam > np.min(w3response.W3)) & (newlam < np.max(w3response.W3))
    lamF_lam_w3 = rebin_spec(np.array(newlam[w3ndx]), \
                          np.array(totalflux[w3ndx]),w3response.W3)
    fluxmw3 = integrate.simps((w3response.RSR*lamF_lam_w3), w3response.W3)#/5.5069

    w4response = pd.read_csv(r'w4response.txt', sep = '\s+')
    w4ndx = (newlam > np.min(w4response.W4)) & (newlam < np.max(w4response.W4))
    lamF_lam_w4 = rebin_spec(np.array(newlam[w4ndx]), \
                          np.array(totalflux[w4ndx]),w4response.W4)
    fluxmw4 = integrate.simps((w4response.RSR*lamF_lam_w4), w4response.W4)#/4.1013

#Use zero-mag normalization based on first column of Table 1, Jarrett+ 11
#https://iopscience.iop.org/article/10.1088/0004-637X/735/2/112/pdf    
    mw1 = -2.5*np.log10(fluxmw1/5.4188e-15)
    mw2 = -2.5*np.log10(fluxmw2/2.51720e-15)
    mw3 = -2.5*np.log10(fluxmw3/3.58781e-16)
    mw4 = -2.5*np.log10(fluxmw4/2.0876e-17)
    
#    print(mw1)
#    print(mw1-mw2)
#    print(mw2-mw3)
    
    w23 = mw2 - mw3
    w12 = mw1 - mw2    
    ax.plot(w23, w12, 'o')
