# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 16:29:16 2019

@author: mugdhapolimera

Reading and plotting SAMI IFU emission line measurements on a BPT plot

"""

import numpy as np
import pandas as pd
from astropy.table import Table as table
from astropy.io import fits
import matplotlib.pyplot as plt
import os

#define demarcation function: log_NII_HA vs. log_OIII_HB
def n2hacompmin(log_NII_HA): #composite minimum line from equation 1, Kewley 2006
    return 1.3 + (0.61 / (log_NII_HA - 0.05))
def n2hamain(log_NII_HA): #main line for NII/H-alpha from equation 5, Kewley 2006
    return 1.19 + (0.61 / (log_NII_HA - 0.47))
def s2hamain(log_SII_HA): #main line for SII/H-alpha from equation 2, Kewley 2006
    return 1.30 + (0.72 / (log_SII_HA - 0.32))
def s2halinseyf(log_SII_HA): #liner/seyfert divider for SII/H-alpha
    return 0.76 + 1.89*log_SII_HA
def o1hamain(log_OI_HA): #main line for OI/H-alpha from equation 3, Kewley 2006
    return 1.33 + (0.73 / (log_OI_HA + 0.59))
def o1halinseyf(log_OI_HA): #liner/seyfert divider for OI/H-alpha
    return 1.3 + 1.18*log_OI_HA
def o1hacrit(log_OI_HA): #boundary for OI/H-alpha
    return -0.59
def ratioerror(num,num_err,den, den_err):
    err = (num/den) * np.sqrt((num_err/num)**2 + (den_err/den)**2)
    return err

#folder = '71146'
#folder = '372320'
folder = '372374'
center = (25,25)
os.chdir(r'C:\Users\mugdhapolimera\Desktop\UNC\Courses\Research\SAMI Data\/'+folder)
#filename = 'adaptive_1-comp.fits'
#filename = 'adaptive_recom-comp.fits'
filename = 'default_recom-comp.fits'
#filename = 'extinct-corr_default_recom-comp.fits'
error = 0
halpha = fits.open(folder+'_Halpha_'+filename)[0].data[0,:,:]#.flatten()
halpha_err = fits.open(folder+'_Halpha_'+filename)[1].data[0,:,:]#.flatten()
hbeta = fits.open(folder+'_Hbeta_'+filename)[0].data#.flatten()
hbeta_err = fits.open(folder+'_Hbeta_'+filename)[1].data#.flatten()
oiii = fits.open(folder+'_OIII5007_'+filename)[0].data#.flatten()
oiii_err = fits.open(folder+'_OIII5007_'+filename)[1].data#.flatten()
oi = fits.open(folder+'_OI6300_'+filename)[0].data#.flatten()
oi_err = fits.open(folder+'_OI6300_'+filename)[0].data#.flatten()
nii = fits.open(folder+'_NII6583_'+filename)[0].data#.flatten()
nii_err = fits.open(folder+'_NII6583_'+filename)[1].data#.flatten()
#sii = (fits.open(folder+'_SII6716_'+filename)[0].data + \
#        fits.open(folder+'_SII6731_'+filename)[0].data**2).flatten()
#sii_err = np.sqrt(fits.open(folder+'_SII6716_'+filename)[0].data**2 + \
#              fits.open(folder+'_SII6731_'+filename)[0].data**2).flatten()
sii_1 = fits.open(folder+'_SII6716_'+filename)[0].data#.flatten()
sii_1_err = fits.open(folder+'_SII6716_'+filename)[0].data#.flatten()
sii_2 = fits.open(folder+'_SII6731_'+filename)[0].data#.flatten()
sii_2_err = fits.open(folder+'_SII6731_'+filename)[0].data#.flatten()

ha_cen = halpha[center]
hb_cen = hbeta[center]
nii_cen = nii[center]
oiii_cen = oiii[center]

good = [~np.isnan(halpha) & ~np.isnan(hbeta) & ~np.isnan(nii) & ~np.isnan(sii_1)
        & ~np.isnan(sii_2) & ~np.isnan(oi) & ~np.isnan(oiii)]
halpha = halpha[good]
hbeta = hbeta[good]
nii = nii[good]
sii_1 = sii_1[good]
sii_2 = sii_2[good]
oi = oi[good]
oiii = oiii[good]

halpha_err = halpha_err[good]
hbeta_err = hbeta_err[good]
nii_err = nii_err[good]
sii_1_err = sii_1_err[good]
sii_2_err = sii_2_err[good]
oi_err = oi_err[good]
oiii_err = oiii_err[good]

catid = [folder]*len(halpha)
data = list(zip(catid, hbeta, hbeta_err,oiii, oiii_err, oi, oi_err,
                halpha, halpha_err, nii, nii_err, sii_1, sii_1_err, sii_2, sii_2_err))
names = ['CATID',
         'h_beta_flux', 'h_beta_flux_err', 
       'oiii_5007_flux', 'oiii_5007_flux_err',
       'oi_6300_flux', 'oi_6300_flux_err', 
       'h_alpha_flux','h_alpha_flux_err',
       'nii_6584_flux', 'nii_6584_flux_err', 
       'sii_6717_flux','sii_6717_flux_err',
       'sii_6731_flux', 'sii_6731_flux_err']

df = pd.DataFrame(data, columns = names)
df.to_csv(folder+'.csv')
sii = sii_1 + sii_2
sii_err = np.sqrt(sii_1_err**2 + sii_2_err**2)
n2ha = nii/halpha
s2ha = sii/halpha
o1ha= oi/halpha
o3hb = oiii/hbeta
n2ha_err = ratioerror(nii, nii_err, halpha, halpha_err)
s2ha_err = ratioerror(sii, sii_err, halpha, halpha_err)
o1ha_err = ratioerror(oi, oi_err, halpha, halpha_err)
o3hb_err = ratioerror(oiii, oiii_err, hbeta, hbeta_err)

refn2ha = np.linspace(-3.0, 0.35)
refoiha = np.linspace(-2.5, -0.4)
refsiiha = np.linspace(-2, 0.3,100)

fig = plt.figure()#'NII Scatter Plot')
ax1 = fig.add_subplot(111)
ax1.set_xlim(-1.5,0.5)
ax1.set_ylim(-1.0,1.0)
ax1.plot(refn2ha, n2hamain(refn2ha), 'k', 
                  label = 'ke01 Theoretical Maximum Starburst Line')
ax1.plot(refn2ha[refn2ha < 0], n2hacompmin(refn2ha[refn2ha < 0]),
                      'k-.', label = 'Ka03 Composite Line')
ax1.plot(n2ha, o3hb, 'k.', alpha = 0.5, markersize = 5)#, label = 'Definite Star Forming')
ax1.plot(nii_cen/ha_cen, oiii_cen/hb_cen, 'ro', alpha = 0.5, markersize = 5)#, label = 'Definite Star Forming')
ax1.set_xlabel(r"$\rm \log([NII]/H\alpha)$", fontsize = 22)
ax1.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
if error: 
    ax1.errorbar(n2ha.flatten(), o3hb.flatten(), xerr = n2ha_err.flatten(),
                            yerr = o3hb_err.flatten(), fmt = 'b.', alpha = 0.5,
                        markersize = 8, mew = 0, label = 'SF-to-AGN', ecolor = 'k')

#SII/OIII plot
fig = plt.figure()#'SII Scatter Plot')
ax2 = fig.add_subplot(111)
ax2.plot(refsiiha, s2hamain(refsiiha), 'k',  label = 'Ke01 Line')
ax2.plot(refsiiha[refsiiha > -0.31], s2halinseyf(refsiiha[refsiiha > -0.31]),
                  'k--', label = 'Liner/Seyfert Division')
ax2.set_xlim(-1.5, 0.5)
ax2.set_ylim(-1.0,1.0)
ax2.set_xlabel(r"$\rm \log([SII]/H\alpha)$", fontsize = 22)
ax2.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
ax2.plot(s2ha, o3hb, 'k.', markersize = 5, 
                    alpha = 0.5, label = 'SF')
if error:
    ax2.errorbar(s2ha.flatten(), o3hb.flatten(), xerr = s2ha_err.flatten(),
                            yerr = o3hb_err.flatten(), fmt = 'b.', alpha = 0.5,
                        markersize = 8, mew = 0, label = 'SF-to-AGN', ecolor = 'k')

#OI/OIII plot
fig = plt.figure()#'OI Scatter Plot')
ax3 = fig.add_subplot(111)
ax3.plot(refoiha[refoiha < -0.7], o1hamain(refoiha[refoiha < -0.7]),
                  'k', label = 'Ke01 Theoretical Maximum Starburst Line')
ax3.plot(refoiha[refoiha < -0.7], o1hamain(refoiha[refoiha < -0.7]),
                  'k-.', label = 'Ka03 Composite Line')
#ax3.set_xlim(-2.0, -0.4)
#ax3.set_ylim(-1.0,1.0)
ax3.plot(refoiha[refoiha > -1.13], o1halinseyf(refoiha[refoiha > -1.13]),
                               'k--', label = 'Ke06 Liner/Seyfert Division Line')
ax3.set_xlabel(r"$\rm \log([OI]/H\alpha)$", fontsize = 22)
ax3.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
ax3.plot(o1ha, o3hb, 'k.', alpha = 0.5, 
                    markersize = 5, label = 'SF')
if error:
    ax3.errorbar(o1ha.flatten(), o3hb.flatten(), xerr = o1ha_err.flatten(),
                                yerr = o3hb_err.flatten(), fmt = 'b.', alpha = 0.5,
                            markersize = 8, mew = 0, label = 'SF-to-AGN', ecolor = 'k')
df = pd.read_csv(folder+'_smcdext.csv')
#create line ratios/H-alpha and [OIII]/H-beta
nii_2 = df['nii_6584_flux']
nii_sum = df['nii_6584_flux']
nii_sum_err = df['nii_6584_flux_err']**2
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
sii_sum = df['sii_6717_flux'] + df['sii_6731_flux']
sii_sum_err = np.sqrt(df['sii_6717_flux_err']**2 + df['sii_6731_flux_err']**2)

#Filter Data: all non-negative SEL fluxes and errors; Hbeta >3sigma
gooddata = ((h_alpha > 0) & (nii_sum > 0) & (oiii > 0) & (oi > 0) &
            (sii_sum > 0) & (h_beta > 0) & (h_beta > 3*h_beta_err) &
            (h_alpha_err > 0) & (nii_sum_err > 0) & (oiii_err > 0) & 
            (oi_err > 0) & (sii_sum_err > 0))

data = gooddata #use ALL galaxy data within catalog
nii = nii[data]
nii_sum = nii_sum[data]
oiii = oiii[data]
oiii_err = oiii_err[data]
oi = oi[data]
oi_err = oi_err[data]
sii_sum = sii_sum[data]
sii_sum_err = sii_sum_err[data]
h_beta = h_beta[data]
h_beta_err = h_beta_err[data]
h_alpha = h_alpha[data]
h_alpha_err = h_alpha_err[data]

fig = plt.figure('NII Scatter Plot - deextincted')
ax1 = fig.add_subplot(111)
ax1.set_xlim(-1.5,0.5)
ax1.set_ylim(-1.0,1.0)
ax1.plot(refn2ha, n2hamain(refn2ha), 'k', 
                  label = 'ke01 Theoretical Maximum Starburst Line')
ax1.plot(refn2ha[refn2ha < 0], n2hacompmin(refn2ha[refn2ha < 0]),
                      'k-.', label = 'Ka03 Composite Line')
ax1.plot(n2ha, o3hb, 'k.', alpha = 0.5, markersize = 5)#, label = 'Definite Star Forming')
ax1.set_xlabel(r"$\rm \log([NII]/H\alpha)$", fontsize = 22)
ax1.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
if error: 
    ax1.errorbar(n2ha.flatten(), o3hb.flatten(), xerr = n2ha_err.flatten(),
                            yerr = o3hb_err.flatten(), fmt = 'b.', alpha = 0.5,
                        markersize = 8, mew = 0, label = 'SF-to-AGN', ecolor = 'k')

#SII/OIII plot
fig = plt.figure('SII Scatter Plot - deextincted')
ax2 = fig.add_subplot(111)
ax2.plot(refsiiha, s2hamain(refsiiha), 'k',  label = 'Ke01 Line')
ax2.plot(refsiiha[refsiiha > -0.31], s2halinseyf(refsiiha[refsiiha > -0.31]),
                  'k--', label = 'Liner/Seyfert Division')
ax2.set_xlim(-1.5, 0.5)
ax2.set_ylim(-1.0,1.0)
ax2.set_xlabel(r"$\rm \log([SII]/H\alpha)$", fontsize = 22)
ax2.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
ax2.plot(s2ha, o3hb, 'k.', markersize = 5, 
                    alpha = 0.5, label = 'SF')
if error:
    ax2.errorbar(s2ha.flatten(), o3hb.flatten(), xerr = s2ha_err.flatten(),
                            yerr = o3hb_err.flatten(), fmt = 'b.', alpha = 0.5,
                        markersize = 8, mew = 0, label = 'SF-to-AGN', ecolor = 'k')

#OI/OIII plot
fig = plt.figure('OI Scatter Plot - deextincted')
ax3 = fig.add_subplot(111)
ax3.plot(refoiha[refoiha < -0.7], o1hamain(refoiha[refoiha < -0.7]),
                  'k', label = 'Ke01 Theoretical Maximum Starburst Line')
ax3.plot(refoiha[refoiha < -0.7], o1hamain(refoiha[refoiha < -0.7]),
                  'k-.', label = 'Ka03 Composite Line')
#ax3.set_xlim(-2.0, -0.4)
#ax3.set_ylim(-1.0,1.0)
ax3.plot(refoiha[refoiha > -1.13], o1halinseyf(refoiha[refoiha > -1.13]),
                               'k--', label = 'Ke06 Liner/Seyfert Division Line')
ax3.set_xlabel(r"$\rm \log([OI]/H\alpha)$", fontsize = 22)
ax3.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
ax3.plot(o1ha, o3hb, 'k.', alpha = 0.5, 
                    markersize = 5, label = 'SF')
if error:
    ax3.errorbar(o1ha.flatten(), o3hb.flatten(), xerr = o1ha_err.flatten(),
                                yerr = o3hb_err.flatten(), fmt = 'b.', alpha = 0.5,
                            markersize = 8, mew = 0, label = 'SF-to-AGN', ecolor = 'k')
