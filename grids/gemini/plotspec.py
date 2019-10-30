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

spectype = 'bpass'
#spectype = 'cloudy'

if spectype == 'bpass':
    spec = pd.read_csv('spectra-bin.z008_1.dat',sep = '\s+', header=None, 
                   usecols = [0,41], names = ['lam', 'flux'])
if spectype == 'cloudy':
    spec = pd.read_csv('agn_0000-Z_0400.con',sep = '\s+', header=None, skiprows = 1,
                   usecols = [0,6], names = ['lam', 'flux'])

#spec = spec.iloc[np.where((spec.lam < 10000 ) & (spec.lam >=1 ))]
spec = spec.iloc[np.where((8000 > spec.lam) & (spec.lam > 4000))]
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
new_spec, vind = pyasl.specAirVacConvert(newlam, new_spec, \
                direction="vactoair")
hbeta_flux = new_spec[np.where(newlam == newlam[newlam>=4864.55][0])]
new_spec = new_spec/hbeta_flux
plt.plot(newlam,new_spec, 'k')

resdata = readsav('../../../SDSS_spectra/resolvecatalog.dat')
galname = 'rf0376'
resphot = readsav('../../../SDSS_spectra/resolvecatalogphot.dat')
galndx = np.where(resdata['name'] == galname)
v = resdata['vhel'][galndx] #km/s
z = v/3e5 #redshift
r_pc = (v/70)*10**6*(1+z) #pc
r = r_pc*3.086e+18 #pc to cm
mstar = resdata['mstars'][np.where(resdata['name'] == galname)]
#new_spec = new_spec/4*np.pi*r
emlines = pd.read_csv('../../../izi/Richardson-0-0_1-0agn-M-5_0-BPASS-Binary-CSF-n=1e3-40.0Myr-NichollsCE-D_G-RR14_Fstar_0_3-unified-2.csv')
emlines = emlines[(emlines.AGNFRAC == 1.0) & (emlines.LOGZ == np.log10(0.4)) &\
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
        sigma = 1.69/2.355 #FWHM to sigma (69/3e5)*6500
        x = (newlam - mu)/sigma
        line = norm.pdf(x)*np.array(emlines[lam])
        #print line, max(line), lam
        line = (line/max(line)) * np.array(emlines[lam]) 
        emspec+= line
        plt.plot(newlam,line,'g')
    if lam == 'halpha':
        #converting broad gaussian width from velocity into wavelength
        sigma = ((500/2.355)/3e5)*mu 
        x = (newlam - mu)/sigma
        line = norm.pdf(x)*0.15*np.array(emlines[lam]) #15% of height of Halpha
        #print line, max(line), lam
        line = (line/max(line))*0.15*np.array(emlines[lam])
        line = line*emlines.AGNFRAC.iloc[0]
        emspec+= line
        plt.plot(newlam,line,'g')

#Fit gaussian to Balmer lines to model the absorption trough
from scipy import optimize

#def gaussian(x, amplitude, stddev):
#    return amplitude * np.exp(-((x - 4861) / 4 / stddev)**2)
def gaussian(x, amplitude, stddev, offset):
    return amplitude * np.exp(-(x - mean) ** 2 / (2 * stddev ** 2)) + offset
mean = wavelengths['hbeta']
hbndx = (newlam > 4800) & (newlam < 4920)
hb = np.array(new_spec[hbndx])
hbetalam = np.array(newlam[hbndx])
popt, _ = optimize.curve_fit(gaussian, hbetalam, hb, p0 = [-100,100,100])
plt.plot(hbetalam,hb,'k.')
plt.plot(hbetalam,gaussian(hbetalam,*popt),'orange')
sig = ((250/2.354)/3e5)*mean #converting 190 km/s to sigma in Angstroms 
sigratio = popt[1]/sig
plt.plot(hbetalam,gaussian(hbetalam,popt[0]*sigratio/4, sig, popt[-1])-popt[-1],'m')
emspec += gaussian(newlam,popt[0]*sigratio/4, sig, popt[-1])-popt[-1]

mean = wavelengths['halpha']
handx = (newlam > 6480) & (newlam < 6630)
ha = np.array(new_spec[handx])
halphalam = np.array(newlam[handx])
popt, _ = optimize.curve_fit(gaussian, halphalam, ha, p0 = [-100,100,100])
plt.plot(halphalam,ha,'k.')
plt.plot(halphalam,gaussian(halphalam,*popt),'orange')
sigratio = popt[1]/sig
plt.plot(halphalam,gaussian(halphalam,popt[0]*sigratio/4, sig, popt[-1])-popt[-1],'m')
#emspec += gaussian(newlam,popt[0]*sigratio/4, sig, popt[-1])-popt[-1]
print(emspec)  
#totalspec = (new_spec/max(new_spec))+(emspec/max(emspec))
lsun = 3.826*10**33 #ergs/s
age = 10**(6+0.1*(42-2)) # in years
#norm =1# 25*1e-2
#totalflux = norm*mstar*(totalspec* lsun) /(4*np.pi*r*r*1e6)
#totalflux = totalspec #new_spec/(newlam)

optndx = (newlam < 7400) & (newlam > 4000)
optlam = newlam[optndx]
optspec = new_spec[optndx] 

notmask = ((optlam < 4070) | (optlam > 4150)) \
        & ((optlam < 4300) | (optlam > 4400)) \
        & ((optlam < 4810) | (optlam > 5070)) \
        & ((optlam < 5830) | (optlam > 5940)) \
        & ((optlam < 6260) | (optlam > 6360)) \
        & ((optlam < 6500) | (optlam > 6800)) \
        & ((optlam < 7000) | (optlam > 7200)) \
        & ((optlam < 7250) | (optlam > 7380))
cont = np.poly1d(np.polyfit(optlam[notmask], optspec[notmask],10))
totalflux = np.zeros(len(new_spec))
totalflux[optndx]  = cont(optlam)
plt.plot(optlam,cont(optlam),'m')
plt.plot(optlam[notmask],optspec[notmask],'go')
fullmask = ((newlam <4070) | (newlam > 4150)) \
        & ((newlam <4300) | (newlam > 4400)) \
        & ((newlam <4810) | (newlam > 5070)) \
        & ((newlam <5830) | (newlam >5940)) \
        & ((newlam < 6260) | (newlam > 6360)) \
        & ((newlam < 6500) | (newlam > 6800)) \
        & ((newlam < 7000) | (newlam > 7200)) \
        & ((newlam < 7250) | (newlam > 7380))
#To match the line/continuum ratio from SDSS for rf0376
eq_width = 15*1e-17#/23.0
totalflux[fullmask] = new_spec[fullmask]
#totalflux = new_spec
totalflux = (totalflux+emspec)*eq_width#hbeta_flux
if spectype == 'cloudy':
    totalflux = (totalflux/(4*np.pi*r*r))/newlam

#if spectype == 'bpass':
    #totalflux = mstar*(totalflux* lsun) /(4*np.pi*r*r*age)

plt.figure()
plt.plot(newlam, totalflux/1e-17)
plt.xlabel(r'Wavelength (in ${\AA}$)')
plt.ylabel(r'Flux Density ($10^{-17}$ $ergs/s/cm^2/{\AA}$)')
plt.xlim(4850,5050)

plt.figure()
plt.plot(newlam, totalflux/1e-17)
plt.xlim(6250,6600)
plt.xlabel(r'Wavelength (in ${\AA}$)')
plt.ylabel(r'Flux Density ($10^{-17}$ $ergs/s/cm^2/{\AA}$)')
#plt.plot(optlam,cont(optlam))
#spectrum = pd.DataFrame(data = {'flux' : totalflux[optndx]*10, 'lam' : newlam[optndx]/10})
#spectrum.to_csv('mock_spectrum.csv')
np.savetxt('mock_spectrum'+str(emlines.AGNFRAC.iloc[0])+'.txt', zip(newlam[optndx]/10,totalflux[optndx]))
#plt.xlim(0,10000)
#add reddening
#Assume E_BV of one of our galaxies?

def func(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))    
def composite_spectrum(x, # data
                       a1, x01, sigma1, # 1st line
                       a2, x02, sigma2, # 2nd line
                       a3, x03, sigma3, # 3rd line
                       a4, x04, sigma4): # 4th line    
    return (func(x, a1, x01, sigma1) \
            + func(x, a2, x02, sigma2) \
            + func(x, a3, x03, sigma3)\
            + func(x, a4, x04, sigma4))
#                    + func(x, a3, x03, sigma3))

guess = [1000, 6563, 1.69/2.354, 500, 6584, 1.69/2.354,
         200, 6548, 1.69/2.354, 150, 6563, 4]

from scipy.optimize import curve_fit
x2 = newlam[(newlam>6500) & (newlam<6600)]
y2 = totalflux[(newlam>6500) & (newlam<6600)]
y2 = y2/np.max(y2)
popt, pcov = curve_fit(composite_spectrum, x2, y2, p0 = guess)
plt.figure()
plt.plot(x2,y2)
plt.plot(x2, composite_spectrum(x2, *popt), 'k', label='Total fit')
plt.plot(x2, func(x2, *popt[-3:]), c='r', 
         label='Broad component')
FWHM = round(2*np.sqrt(2*np.log(2))*popt[-1],4)
#plt.axvspan(popt[-2]-FWHM/2, popt[-2]+FWHM/2, 
#            facecolor='g', alpha=0.3, label='FWHM = %s'%(FWHM))
plt.legend(fontsize=10)
#plt.ylim(0,1.1)
plt.show()

import numpy as np
from scipy.optimize import leastsq
import matplotlib.pyplot as plt

######################################
# Setting up test data
def norm(x, mean, sd):
  norm = []
  for i in range(x.size):
    norm += [1.0/(sd*np.sqrt(2*np.pi))*np.exp(-(x[i] - mean)**2/(2*sd**2))]
  return np.array(norm)

mean1, mean2 = 0, -2
std1, std2 = 0.5, 1 

x = x2 #np.linspace(-20, 20, 500)
y_real = y2# norm(x, mean1, std1) + norm(x, mean2, std2)

######################################
# Solving
m, dm, sd1, sd2 = [6563, 0, 2, 10]
p = [m, dm, sd1, sd2] # Initial guesses for leastsq
y_init = norm(x, m, sd1) + norm(x, m+dm, sd2) # For final comparison plot

def res(p, y, x):
  m, dm, sd1, sd2 = p
  m1 = m
  m2 = m1 + dm
  y_fit = norm(x, m1, sd1) + norm(x, m2, sd2)
  err = y - y_fit
  return err

plsq = leastsq(res, p, args = (y_real, x))

y_est = norm(x, plsq[0][0], plsq[0][2]) + norm(x, plsq[0][0] + plsq[0][1], plsq[0][3])
plt.figure()
plt.plot(x, y_real, label='Real Data')
plt.plot(x, y_init, 'r.', label='Starting Guess')
plt.plot(x, y_est, 'g', label='Fitted')
plt.legend()
plt.show()

fwhminpix=[]
arcspec = y2
lines = [6548, 6563, 6561, 6584]
waverange = x2
add=[10, 10, 40, 10]

shap=np.shape(arcspec)

conv=((waverange[1]-waverange[0])/(shap[0]-1))

lines=np.array(lines)
pixlines=(lines/conv)-(waverange[0]/conv)
pixlines = [int(i) for i in pixlines]
from mpfit import mpfit

for i in np.arange(len(lines)):
    nput=arcspec[shap[0]//2+100,pixlines[i]-add[i]:pixlines[i]+add[i]]
    gerr=np.zeros(len(nput))+0.5
    xaxis=np.arange(len(nput))

    max=np.max(nput)
    mid=np.argmax(nput)
    print(i,mid,max)

    p0=[max,mid,5.,30.]

    def myfunct(p, fjac=None, x=None, y=None, err=None):
        model = p[0] * np.exp(-((x-p[1])**2.)/(2.*p[2]**2.)) + p[3]
        status = 0
        return([status, (y-model)/err])

    fa = {'x':xaxis, 'y':nput, 'err':gerr}

    m=mpfit(myfunct,p0,functkw=fa)

    def func2(x, a, b, d, c):
        return a * np.exp(-((x-b)**2.)/(2.*d**2.)) + c

    fitdata=func2(xaxis,m.params[0],m.params[1],m.params[2],m.params[3])

    #figg1=ax7.add_subplot(3,1,i+1)
    sigma=m.params[2]
    fwhm=np.abs(sigma*2.35482)
    fwhminpix.append(fwhm)
    
    ax5=plt.subplot2grid((2,5),(1,i))#,colspan=4)
    plt.plot(xaxis,nput)
    plt.plot(xaxis,fitdata)
    plt.title('gaussian fits to arc lines')
    plt.ylim([np.min(nput),np.max(nput)])
    
    ax5.annotate('FWHM =',xy=(np.min(xaxis)+3.,np.max(nput)/2.0))
    ax5.annotate('%3.3g' % (fwhm*conv) + ' A',xy=(np.min(xaxis)+4.,np.max(nput)/2.3))

fwhminpix=np.array(fwhminpix)

fwhminA=fwhminpix*conv
