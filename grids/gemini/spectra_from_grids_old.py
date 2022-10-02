# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 21:53:38 2022

@author: mugdhapolimera
"""

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

#spectype = 'bpass'
spectype = 'cloudy'

folder = r'C:\Users\mugdhapolimera\github\SDSS_spectra\mid_ir\satyapalmodels\/'
if spectype == 'bpass':
    spec = pd.read_csv('spectra-bin.z008_1.dat',sep = '\s+', header=None, 
                   usecols = [0,41], names = ['lam', 'flux'])
if spectype == 'cloudy':
#    speclowmbh = pd.read_csv('agn_0320-Z_0200_M-3_0.con',sep = '\s+', header=None, skiprows = 1,
#                       usecols = [0,6], names = ['lam', 'flux']) #agn_0000-Z_0400.con
#    spechighmbh = pd.read_csv('agn_0320-Z_1000_M-8_0.con',sep = '\s+', header=None, skiprows = 1,
#                   usecols = [0,6], names = ['lam', 'flux']) #agn_0000-Z_0400.con

    speclowmbh = pd.read_csv(folder+r'grid000000002_Z_1_n_300_100per.con',sep = '\s+', header=None, skiprows = -2,
                       usecols = [0,6], names = ['lam', 'flux'])[:-1] #agn_0000-Z_0400.con
    speclowmbh = speclowmbh.astype(np.float64)
    spechighmbh = pd.read_csv(folder+r'grid000000013_Z_1_n_300_100per.con',sep = '\s+', header=None, skiprows = 1,
                   usecols = [0,6], names = ['lam', 'flux'])[:-1] #agn_0000-Z_0400.con
    spechighmbh = spechighmbh.astype(np.float64)


bhspecfiles = glob.glob(folder+'*.con')


#specs = [speclowmbh, spechighmbh]

figspec, specax = plt.subplots(1)
specax.set_xlabel(r'Wavelength (in ${\mu m}$)')
specax.set_ylabel(r'Flux Density ')#($10^{-17}$ $ergs/s/cm^2/{\AA}$)')
specax.set_yscale('log')
specax.set_xscale('log')
#np.savetxt('mock_spectrum'+str(emlines.AGNFRAC.iloc[0])+'.txt', zip(newlam[optndx]/10,totalflux[optndx]*1e35))
specax.set_xlim(0.1,100)
#    specax.ylim(0.001,10)
specax.axvspan(3,3.8, alpha = 0.3, color = 'k')
specax.axvspan(4,5, alpha = 0.3, color = 'k')
specax.axvspan(9,13, alpha = 0.3, color = 'k')
specax.axvspan(20,25, alpha = 0.3, color = 'k')

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

#for spec in specs:
for specfile in bhspecfiles:
    spec = pd.read_csv(specfile,sep = '\s+', header=None, skiprows = -2,
                       usecols = [0,6], names = ['lam', 'flux'])[:-1] 
    spec = spec.astype(np.float64)
    spec = spec.iloc[np.where((spec.lam < 1e6 ) & (spec.lam >=1000 ))]
    #spec = spec.iloc[np.where((8000 > spec.lam) & (spec.lam > 4000))]
    #plt.figure()
    #plt.plot(spec.lam,spec.flux/np.max(spec.flux))
    newlam = np.arange(np.min(spec.lam),np.max(spec.lam),0.5)
    spec.flux = np.array(spec.flux[::-1])
    spec.lam = np.array(spec.lam[::-1])
    new_spec = rebin_spec(np.array(spec.lam),np.array(spec.flux),newlam)
    #plt.plot(newlam,new_spec/np.max(new_spec), 'r')
    
    #vacuum to air conversion
    s = 10**4/newlam
    n = 1 + 0.0000834254 + 0.02406147 / (130 - s**2) + 0.00015998 / (38.9 - s**2)
    airlam = newlam/n 
    #air_spec = rebin_spec(np.array(newlam),np.array(new_spec),airlam)
    new_spec = new_spec*n
    newlam = airlam
    
    #plt.plot(airlam,new_spec/max(new_spec), 'k')
    #new_spec, vind = pyasl.specAirVacConvert(newlam, new_spec, direction="vactoair")
    hbeta_flux = new_spec[np.where(newlam == newlam[newlam>=4864.55][0])]
    new_spec = new_spec#/hbeta_flux
#    plt.plot(newlam,new_spec, 'k')
    
    resdata = readsav('../../../SDSS_spectra/resolvecatalog.dat')
    galname = 'rs0010'
    resphot = readsav('../../../SDSS_spectra/resolvecatalogphot.dat')
    galndx = np.where(resdata['name'] == galname)
    v = resdata['vhel'][galndx] #km/s
    z = v/3e5 #redshift
    r_pc = (v/70)*1e6 *(1+z) #pc
    r = np.float64(r_pc*3.086e+18) #pc to cm
    #mstar = resdata['mstars'][np.where(resdata['name'] == galname)]
    ##new_spec = new_spec/4*np.pi*r
    #sdss = 0
    #if sdss:
    #    #emlines = pd.read_csv('../../../SDSS_spectra/RESOLVE_full_blend_dext_new.csv')
    #    emlines = pd.read_csv('../../../SDSS_spectra/RESOLVE_full_snr5.csv')
    #    emlines.index = emlines.name
    #    emlines = emlines.loc['rs0010']
    #    wavelengths = {'oii_3726_flux' : 3726.032, 'oii_3729_flux' : 3728.815, 
    #                   'neiii_3869_flux' : 3868.760, 'h_gamma_flux' : 4340.471, 
    #                   'oiii_4363_flux' : 4363.210, 'heii_4685_flux' : 4685.710,
    #                   'h_beta_flux' : 4861.333, 'oiii_4959_flux': 4958.911, 
    #                   'oiii_5007_flux' : 5006.843, 'hei_5875_flux' : 5875.624, 
    #                   'oi_6300_flux' : 6300.304, 'nii_6548_flux' : 6548.050, 
    #                   'h_alpha_flux' : 6562.819, 'nii_6584_flux' : 6583.46, 
    #                   'sii_6717_flux' : 6716.440, 'sii_6731_flux' : 6730.810,
    #                   'ariii_7136_flux' : 7135.790}
    #else:
    #    emlines = pd.read_csv('../../../izi/Richardson-0-0_1-0agn-M-5_0-BPASS-Binary-CSF-n=1e3-40.0Myr-NichollsCE-D_G-RR14_Fstar_0_3-unified.csv')
    #    emlines = emlines[(emlines.AGNFRAC == 0.5) & (emlines.LOGZ == np.log10(0.4)) &\
    #                       (7.3>emlines.LOGQ) & (emlines.LOGQ > 7.2)]
    #    wavelengths = {'oii3726' : 3726.032, 'oii3729' : 3728.815, 
    #                   'neiii3869' : 3868.760, 'hgamma' : 4340.471, 
    #                   'oiii4363' : 4363.210, 'heii4685' : 4685.710,
    #                   'hbeta' : 4861.333, 'oiii4959': 4958.911, 
    #                   'oiii5007' : 5006.843, 'hei5875' : 5875.624, 
    #                   'oi6300' : 6300.304, 'nii6548' : 6548.050, 
    #                   'halpha' : 6562.819, 'nii6584' : 6583.46, 
    #                   'sii6717' : 6716.440, 'sii6731' : 6730.810,
    #                   'ariii7136' : 7135.790}
    #emspec = np.zeros(len(new_spec))
    #sdss_hbeta = emlines['hbeta']#_flux']
    ##plt.figure()
    #for lam in wavelengths.keys():
    #    mu = wavelengths[lam]
    #    if mu > newlam[0]:
    #        if lam in emlines.keys():
    #            if sdss:
    #                emlines[lam] = emlines[lam]/sdss_hbeta
    #            sigma = 1.69/2.355 #FWHM to sigma (69/3e5)*6500
    #            x = (newlam - mu)/sigma
    #            line = norm.pdf(x)*np.array(emlines[lam])
    #            #print line, max(line), lam
    #            line = (line/np.max(line)) * np.array(emlines[lam]) 
    #            emspec+= line
    #            plt.plot(newlam,line,'g')
    #    if (lam == 'halpha'):
    #        #converting broad gaussian width from velocity into wavelength
    #        sigma = ((500/2.355)/3e5)*mu 
    #        x = (newlam - mu)/sigma
    #        line = norm.pdf(x)*0.15*np.array(emlines[lam]) #15% of height of Halpha
    #        #print line, max(line), lam
    #        line = (line/np.max(line))*0.15*np.array(emlines[lam])
    #        line = line*emlines.AGNFRAC.iloc[0]
    #        emspec+= line
    #        plt.plot(newlam,line,'g')
    #
    #eq_width = 15*1e-17#/23.0
    #totalflux = (new_spec+emspec)#*eq_width#hbeta_flux
    ##lsun = 3.826*10**33 #ergs/s
    ##age = 10**(6+0.1*(42-2)) # in years
    ##
    ##if spectype == 'cloudy':
    ##    totalflux = (totalflux/(4*np.pi*r*r))/newlam
    ##
    ##if spectype == 'bpass':
    ##    totalflux = mstar*(totalflux* lsun) /(4*np.pi*r*r*age)
    totalflux = new_spec#/(4*np.pi*r*r) #units of ergs/s/cm^2 -- nu*F_nu
#    totalflux = totalflux/(3e8 / newlam*1e-6) #convert nu_fnu to Fnu
    totalflux = totalflux*1e-7 #convert to W/cm^2

    newlam = newlam*1e-10/1e-6 #convert angstrom to microns
    totalflux = totalflux/newlam #CONVERT LAM*f_LAM to f_lam
    #to redshift the spectrum
#    totalflux = totalflux/(1+z) #lam = lam*(1+z); flux = flux/(1+z); flux density lam*F_lam stays same
#    newlam = newlam*(1+z)
    

    
#    nu = 3e8 / newlam*1e-6 #nu = c/lambda
#    totalflux = totalflux/nu #convert nu_fnu to Fnu    
#    totalflux = totalflux / 1e-23 #convert F_nu to Jansky
    specax.plot(newlam, totalflux) #/1e-17)
    
    from scipy import integrate
    w1ndx = (newlam > 3) & (newlam < 3.8)
    w2ndx = (newlam > 4) & (newlam < 5)
    w3ndx = (newlam > 9) & (newlam < 13)
    w4ndx = (newlam > 20) & (newlam < 25)

    w1response = pd.read_csv(r'C:\Users\mugdhapolimera\github\SDSS_spectra\mid_ir\w1response.txt', sep = '\s+')
    w1ndx = (newlam > np.min(w1response.W1)) & (newlam < np.max(w1response.W1))
    lamF_lam_w1 = rebin_spec(np.array(newlam[w1ndx]), \
                          np.array(totalflux[w1ndx]),w1response.W1)
    fluxmw1 = integrate.simps((w1response.RSR*lamF_lam_w1), w1response.W1)#/6.6e-1

    w2response = pd.read_csv(r'C:\Users\mugdhapolimera\github\SDSS_spectra\mid_ir\w2response.txt', sep = '\s+')
    w2ndx = (newlam > np.min(w2response.W2)) & (newlam < np.max(w2response.W2))
    lamF_lam_w2 = rebin_spec(np.array(newlam[w2ndx]), \
                          np.array(totalflux[w2ndx]),w2response.W2)
    fluxmw2 = integrate.simps((w2response.RSR*lamF_lam_w2), w2response.W2)#/1.0423

    w3response = pd.read_csv(r'C:\Users\mugdhapolimera\github\SDSS_spectra\mid_ir\w3response.txt', sep = '\s+')
    w3ndx = (newlam > np.min(w3response.W3)) & (newlam < np.max(w3response.W3))
    lamF_lam_w3 = rebin_spec(np.array(newlam[w3ndx]), \
                          np.array(totalflux[w3ndx]),w3response.W3)
    fluxmw3 = integrate.simps((w3response.RSR*lamF_lam_w3), w3response.W3)#/5.5069

    w4response = pd.read_csv(r'C:\Users\mugdhapolimera\github\SDSS_spectra\mid_ir\w4response.txt', sep = '\s+')
    w4ndx = (newlam > np.min(w4response.W4)) & (newlam < np.max(w4response.W4))
    lamF_lam_w4 = rebin_spec(np.array(newlam[w4ndx]), \
                          np.array(totalflux[w4ndx]),w4response.W4)
    fluxmw4 = integrate.simps((w4response.RSR*lamF_lam_w4), w4response.W4)#/4.1013

    #mw1 = 20.7-2.5*np.log10(fluxmw1)
    #mw2 = 19.5-2.5*np.log10(fluxmw2)
    #mw3 = 17.8-2.5*np.log10(fluxmw3)
    #mw4 = 12.9-2.5*np.log10(fluxmw4)
    #Use Fnu zero points from 
    #https://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html
#    fluxmw1_jy = fluxmw1/1e-23    
#    mw1 = -2.5*np.log10(fluxmw1/309.54)+0.03
#    mw2 = -2.5*np.log10(fluxmw2/171.787)+0.03
#    mw3 = -2.5*np.log10(fluxmw3/31.674)+0.03
#    mw4 = -2.5*np.log10(fluxmw4/8.363)+0.03

#    mw1 = -2.5*np.log10(fluxmw1/8.1787e-15)+0.03
#    mw2 = -2.5*np.log10(fluxmw2/2.4150e-15)+0.03
#    mw3 = -2.5*np.log10(fluxmw3/6.5151e-17)+0.03
#    mw4 = -2.5*np.log10(fluxmw4/5.0901e-18)+0.03

    mw1 = -2.5*np.log10(fluxmw1/5.4188e-15)#+0.03
    mw2 = -2.5*np.log10(fluxmw2/2.51720e-15)#+0.03
    mw3 = -2.5*np.log10(fluxmw3/3.58781e-16)#+0.03
    mw4 = -2.5*np.log10(fluxmw4/2.0876e-17)#+0.03
    
    print(mw1)
    print(mw1-mw2)
    print(mw2-mw3)
    
    w23 = mw2 - mw3
    w12 = mw1 - mw2    
    ax.plot(w23, w12, 'o')
