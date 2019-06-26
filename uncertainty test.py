# -*- coding: utf-8 -*-
"""
Created on Fri May 31 10:52:29 2019

@author: mugdhapolimera
"""

import uncertainties
from uncertainties import unumpy
import pandas as pd
import os 
import numpy as np
import matplotlib.pyplot as plt

os.chdir('C:\Users\mugdhapolimera\github\SDSS_spectra')
resolve = pd.read_csv('RESOLVE_filter_new.csv')

ha = unumpy.uarray(resolve.h_alpha_flux, resolve.h_alpha_flux_err)
oi = unumpy.uarray(resolve.oi_6300_flux, resolve.oi_6300_flux_err)
sii = unumpy.uarray(resolve.sii_6717_flux+resolve.sii_6731_flux,
    np.sqrt(resolve.sii_6717_flux_err**2 + resolve.sii_6731_flux_err**2))

oiii = unumpy.uarray(resolve.oiii_5007_flux, resolve.oiii_5007_flux_err)
hb = unumpy.uarray(resolve.h_beta_flux, resolve.h_beta_flux_err)

oiha = unumpy.log10(oi/ha)
siiha = unumpy.log10(sii/ha)
oiiihb = unumpy.log10(oiii/hb)

oiha_errs = unumpy.std_devs(oiha)
siiha_errs = unumpy.std_devs(siiha)
oiiihb_errs = unumpy.std_devs(oiiihb)

plt.figure()
plt.hist(oiha_errs, bins = 'fd', histtype = 'step', color = 'b')
plt.hist(oiiihb_errs, bins = 'fd', histtype = 'step', color = 'r')
plt.hist(siiha_errs, bins = 'fd', histtype = 'step', color = 'g')

print np.median(oiha_errs),np.median(siiha_errs),np.median(oiiihb_errs)