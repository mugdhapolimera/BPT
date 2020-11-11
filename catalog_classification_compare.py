# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 13:42:37 2020

@author: mugdhapolimera

Compare classifications NSA catalog to JHU and Portsmouth
"""

import pandas as pd
import os
import numpy as np
os.chdir('C:\Users\mugdhapolimera\github\SDSS_spectra')
jhuold = pd.read_csv("RESOLVE_full_snr5_jhu_new.csv")

jhu= pd.read_csv("RESOLVE_full_snr5_dext_jhu.csv")
jhu.index = jhu.name
df = pd.read_csv("RESOLVE_full_snr5_jhu.csv")
df.index = df.name
diff = list(set(jhu.name) ^ set(jhuold.name))

nii = df['nii_6584_flux']
if 'nii_6548_flux' in df.keys():
    nii_sum = (df['nii_6584_flux']+ df['nii_6548_flux'])*3./4
    nii_sum_err = (np.sqrt(df['nii_6584_flux_err']**2 + df['nii_6548_flux_err']**2))*3./4
else:
    nii_sum = df['nii_6584_flux']
    nii_sum_err = df['nii_6584_flux_err']
# note that the ratio uses only the stronger line, but for S/N reasons we add
# the weaker and multiply by 3/4 since Chris Richardson says the canonical
# line ratio is 3:1 (this needs to be updated with a more precise number)
oiii = df['oiii_5007_flux']
oiii_err = df['oiii_5007_flux_err']
h_alpha = df['h_alpha_flux']
h_alpha_err = df['h_alpha_flux_err']
h_beta = df['h_beta_flux']
h_beta_err = df['h_beta_flux_err']
oi = df['oi_6300_flux']
oi_err = df['oi_6300_flux_err']
if 'sii_6717_flux' in df.keys():
    sii_sum = df['sii_6717_flux'] + df['sii_6731_flux']

    sii_sum_err = np.sqrt(df['sii_6717_flux_err']**2 + df['sii_6731_flux_err']**2)
else:
    sii_sum = df['sii_6731_flux']

    sii_sum_err = df['sii_6731_flux_err']
floor = 10**-3
ceil = 1e5
cut = 5
snr = ((h_alpha > cut*h_alpha_err) 
            & (nii_sum > cut*nii_sum_err) & 
           (oiii > cut*oiii_err) & (oi > cut*oi_err) & (sii_sum > cut*sii_sum_err) 
           & (h_beta > cut*h_beta_err))
gooddata = ((h_alpha > floor) & (nii_sum > floor) & (oiii > floor) & (oi > floor) &
                (sii_sum > floor) & (h_beta > floor)  & (h_alpha_err > floor) & 
                (nii_sum_err > floor) & (oiii_err > floor) & 
                (oi_err > floor) & (sii_sum_err > floor) & 
                (h_alpha < ceil) & (nii_sum < ceil) & (oiii < ceil) & (oi < ceil) &
                (sii_sum < ceil) & (h_beta < ceil) & 
                ~np.isnan(h_alpha) & ~np.isnan(nii_sum) & ~np.isnan(oiii) & 
                ~np.isnan(oi) & ~np.isnan(sii_sum) & ~np.isnan(h_beta))
for gal in diff:
    print gal, snr.loc[gal], oi.loc[gal]/oi_err.loc[gal], \
    jhu.oi_6300_flux.loc[gal]/jhu.oi_6300_flux_err.loc[gal]
    