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
from matplotlib import colors
import os
from scipy.io import readsav

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
    #err = (num/den) * np.sqrt((num_err/num)**2 + (den_err/den)**2)
    #err = ((num_err/num) + (den_err/den))/np.log(10)
    err_num2 = (num_err/(num*np.log(10)))**2
    err_den2 = (den_err/(den*np.log(10)))**2
    #err = np.sqrt((num_err/(num*np.log(10)))**2 + (den_err/(den*np.log(10)))**2)
    return np.sqrt(err_num2 + err_den2)
def logerr(err,number):
    log_err = 0.434*(err/number)
    return log_err

resdata = readsav('../SDSS_spectra/resolvecatalog.dat')
galname = 'rs0010'
#galname = 'rs0775'
galname = 'rs0756'
if galname == 'rs0775':
    folder = '71146' #rs0775
if galname == 'rs0010':
    folder = '372320' #rs0010
if galname == 'rs0013':
    folder = '372374'
if galname == 'rs0756':
    folder = '85507'
#folder = '372374' #rs0013 not SNR 5
os.chdir(r'C:\Users\mugdhapolimera\Desktop\UNC\Courses\Research\SAMI Data\/'+ folder)
#filename = 'adaptive_1-comp.fits'
#filename = 'adaptive_recom-comp.fits'
filename = 'default_recom-comp.fits'
#filename = 'extinct-corr_default_recom-comp.fits'
error = 1
halpha = fits.open(folder+'_Halpha_'+filename)[0].data[0,:,:]#.flatten()
halpha_err = fits.open(folder+'_Halpha_'+filename)[1].data[0,:,:]#.flatten()
hbeta = fits.open(folder+'_Hbeta_'+filename)[0].data#.flatten()
hbeta_err = fits.open(folder+'_Hbeta_'+filename)[1].data#.flatten()
oiii = fits.open(folder+'_OIII5007_'+filename)[0].data#.flatten()
oiii_err = fits.open(folder+'_OIII5007_'+filename)[1].data#.flatten()
oi = fits.open(folder+'_OI6300_'+filename)[0].data#.flatten()
oi_err = fits.open(folder+'_OI6300_'+filename)[1].data#.flatten()
nii = fits.open(folder+'_NII6583_'+filename)[0].data#.flatten()
nii_err = fits.open(folder+'_NII6583_'+filename)[1].data#.flatten()
#sii = (fits.open(folder+'_SII6716_'+filename)[0].data + \
#        fits.open(folder+'_SII6731_'+filename)[0].data**2).flatten()
#sii_err = np.sqrt(fits.open(folder+'_SII6716_'+filename)[0].data**2 + \
#              fits.open(folder+'_SII6731_'+filename)[0].data**2).flatten()
sii_1 = fits.open(folder+'_SII6716_'+filename)[0].data#.flatten()
sii_1_err = fits.open(folder+'_SII6716_'+filename)[1].data#.flatten()
sii_2 = fits.open(folder+'_SII6731_'+filename)[0].data#.flatten()
sii_2_err = fits.open(folder+'_SII6731_'+filename)[1].data#.flatten()
sii = sii_1 + sii_2
#pixelscale = 0.5 #arcsec/pix
cube_hdu = fits.open(folder+'_cube_red.fits')[0].header
bluecube_hdu = fits.open(folder+'_cube_blue.fits')[0].header
cube = fits.open(folder+'_cube_red.fits')[0].data#.flatten()
bluecube = fits.open(folder+'_cube_blue.fits')[0].data#.flatten()
pixelscale = cube_hdu['CDELT1'] #deg/pix
v = resdata['vhel'][np.where(resdata['name'] == galname)] #km/s
z = v/3e5
d = (v/70)*10**6 #Hubble's constant
#pc = 2*3.14*d*pixelscale/360 #pc/pix
pc = pixelscale *3600 #in arcsec
MAX = np.max(halpha[~np.isnan(halpha)])
center = [[25],[25]]#np.where(halpha == MAX)
#center = [center[0][0], center[1][0]]
dx = np.arange(np.shape(halpha)[0])- center[0]
dy = np.arange(np.shape(halpha)[1])- center[1]
dxy = np.array(zip(dx,dy))
r = np.zeros(halpha.shape)
for i in range(len(dx)):
    for j in range(len(dy)):
        r[i][j] = np.sqrt((pc*dx[i])**2 + (pc*-dy[j])**2)
rndx = np.where(r<=2.05)
sdss = cube[67:1990,rndx[0],rndx[1]]
nans = np.where(np.isnan(sdss))
sdss[nans] = 0
sdss2arc = np.sum(sdss,axis=1)

hdu = cube_hdu
lam0 = cube_hdu['CRVAL3']-((cube_hdu['CRPIX3']-1)*cube_hdu['CDELT3'])
lam = (np.arange(cube_hdu['NAXIS3']))*cube_hdu['CDELT3'] + lam0
bluelam0 = bluecube_hdu['CRVAL3']-((bluecube_hdu['CRPIX3']-1)*bluecube_hdu['CDELT3'])
bluelam = (np.arange(bluecube_hdu['NAXIS3']))*bluecube_hdu['CDELT3'] + bluelam0
plt.figure()
plt.plot(lam,cube[:,25,25])
plt.plot(bluelam,bluecube[:,25,25])

sdss2arc = sdss2arc.T[:,None,None]
hdu['NAXIS'] = 3
hdu['NAXIS1'] = sdss2arc.shape[1]
hdu['NAXIS2'] = sdss2arc.shape[2]
hdu['NAXIS3'] = sdss2arc.shape[0]

#hdu['CRVAL3'] = hdu['CRVAL1']
hdu['CRVAL1'] = 1
hdu['CRVAL2'] = 1

hdu['CRPIX3'] = 1
hdu['CRVAL3'] = lam[67]

#hdu['CDELT3'] = hdu['CDELT1']
#hdu['CD3_3'] = hdu['CD1_1']
hdu['CDELT1'] = 1
hdu['CDELT2'] = 1


mesh = np.meshgrid(np.arange(50), np.arange(50))
#To search for a particular wavelength
ndx = np.where(abs(lam-6685) == min(abs(lam-6685)))[0][0]
fig, ax = plt.subplots(1)
ax.imshow(cube[ndx,:,:],norm = colors.Normalize(vmin = 0, vmax = 0.04439))
cmap = colors.ListedColormap(['purple','orange','blue','yellow','green','red',\
                              'cyan','magenta','magenta','magenta','magenta',\
                              'magenta','magenta','magenta','magenta','magenta'])
#plt.get_cmap('rainbow',7)#int(np.max(r)))
cax = ax.scatter(mesh[0],mesh[1], s = 3, c = r, cmap = cmap)
fig.colorbar(cax,extend = 'min')

ha_cen = halpha[center]
hb_cen = hbeta[center]
nii_cen = nii[center]
sii_cen = sii[center]
oi_cen = oi[center]
oiii_cen = oiii[center]
ha_cen_err = halpha_err[center]
hb_cen_err = hbeta_err[center]
nii_cen_err = nii_err[center]
#sii_cen_err = sii_err[center]
oi_cen_err = oi_err[center]
oiii_cen_err = oiii_err[center]
center = [center[0][0],center[1][0]]
#halpha = halpha[center[0]-10:center[0]+10,center[1]-10:center[1]+10]
#hbeta = hbeta[center[0]-10:center[0]+10,center[1]-10:center[1]+10]
#nii   = nii[center[0]-10:center[0]+10,center[1]-10:center[1]+10]
#oi    = oi[center[0]-10:center[0]+10,center[1]-10:center[1]+10]
#oiii  = oiii[center[0]-10:center[0]+10,center[1]-10:center[1]+10]
#sii_1   = sii_1[center[0]-10:center[0]+10,center[1]-10:center[1]+10]
#sii_2   = sii_2[center[0]-10:center[0]+10,center[1]-10:center[1]+10]
#halpha_err = halpha_err[center[0]-10:center[0]+10,center[1]-10:center[1]+10]
#hbeta_err = hbeta_err[center[0]-10:center[0]+10,center[1]-10:center[1]+10]
#nii_err   = nii_err[center[0]-10:center[0]+10,center[1]-10:center[1]+10]
#oi_err    = oi_err[center[0]-10:center[0]+10,center[1]-10:center[1]+10]
#oiii_err  = oiii_err[center[0]-10:center[0]+10,center[1]-10:center[1]+10]
#sii_1_err   = sii_1_err[center[0]-10:center[0]+10,center[1]-10:center[1]+10]
#sii_2_err   = sii_2_err[center[0]-10:center[0]+10,center[1]-10:center[1]+10]
#r = r[center[0]-10:center[0]+10,center[1]-10:center[1]+10]
#oi_full = oi
#halpha_full = halpha
#center = [[2],[4]]#np.where(halpha == MAX)
##center = [center[0][0], center[1][0]]
#ha_cen = halpha[center]
#hb_cen = hbeta[center]
#nii_cen = nii[center]
#sii_cen = sii_1[center]+sii_2[center]
#oi_cen = oi[center]
#oiii_cen = oiii[center]

nans = [~np.isnan(halpha) & ~np.isnan(hbeta) & ~np.isnan(nii) & ~np.isnan(sii_1)
        & ~np.isnan(sii_2) & ~np.isnan(oi) & ~np.isnan(oiii)]
err = [(nii/nii_err > 5) & (sii_1/sii_1_err > 5) & (sii_2/sii_2_err > 5) & 
       (oi/oi_err > 5) & (oiii/oiii_err > 5) & (hbeta/hbeta_err > 5) & 
       (halpha/halpha_err > 5)]
sdss = (r <= 2.05)
#good = (nans and err and sdss)
good = pix
halpha = halpha[good]
hbeta = hbeta[good]
nii = nii[good]
sii_1 = sii_1[good]
sii_2 = sii_2[good]
oi = oi[good]
oiii = oiii[good]
r = r[good]
halpha_err = halpha_err[good]
hbeta_err = hbeta_err[good]
nii_err = nii_err[good]
sii_1_err = sii_1_err[good]
sii_2_err = sii_2_err[good]
oi_err = oi_err[good]
oiii_err = oiii_err[good]

good = (r <= 5.0)
halpha = halpha[good]
hbeta = hbeta[good]
nii = nii[good]
sii_1 = sii_1[good]
sii_2 = sii_2[good]
oi = oi[good]
oiii = oiii[good]
r = r[good]
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
n2ha = np.log10(nii/halpha)
s2ha = np.log10(sii/halpha)
o1ha= np.log10(oi/halpha)
o3hb = np.log10(oiii/hbeta)
#n2ha_err = logerr(ratioerror(nii, nii_err, halpha, halpha_err),n2ha)
#s2ha_err = logerr(ratioerror(sii, sii_err, halpha, halpha_err),s2ha)
#o1ha_err = logerr(ratioerror(oi, oi_err, halpha, halpha_err),o1ha)
#o3hb_err = logerr(ratioerror(oiii, oiii_err, hbeta, hbeta_err),o3hb)
n2ha_err = ratioerror(nii, nii_err, halpha, halpha_err)
s2ha_err = ratioerror(sii, sii_err, halpha, halpha_err)
o1ha_err = ratioerror(oi, oi_err, halpha, halpha_err)
o3hb_err = ratioerror(oiii, oiii_err, hbeta, hbeta_err)

refn2ha = np.linspace(-3.0, 0.35)
refoiha = np.linspace(-2.5, -0.4)
refsiiha = np.linspace(-2, 0.32,100)
xlims = [-1.0,0.0]
ylims = [-0.2,0.35]
cmap = plt.get_cmap('rainbow',int(np.max(r)/0.5))
fig,(ax1,ax2,ax3) = plt.subplots(1,3,sharey = True)#'NII Scatter Plot')
ax1.plot(refn2ha, n2hamain(refn2ha), 'k', 
                  label = 'ke01 Theoretical Maximum Starburst Line')
ax1.plot(refn2ha[refn2ha < 0], n2hacompmin(refn2ha[refn2ha < 0]),
                      'k-.', label = 'Ka03 Composite Line')
#ax1.plot(n2ha, o3hb, 'k.', alpha = 0.5, markersize = 5)#, label = 'Definite Star Forming')
ax1.set_xlabel(r"$\rm \log([NII]/H\alpha)$", fontsize = 22)
ax1.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)

#SII/OIII plot
#fig,ax2 = plt.subplots(322)#'SII Scatter Plot')
ax2.plot(refsiiha, s2hamain(refsiiha), 'k',  label = 'Ke01 Line')
ax2.plot(refsiiha[refsiiha > -0.32], s2halinseyf(refsiiha[refsiiha > -0.32]),
                  'k--', label = 'Liner/Seyfert Division')
#ax2.plot(s2ha, o3hb, 'k.', markersize = 5, \
#                    alpha = 0.5, label = 'SF')
#ax2.plot(np.log10(sii_cen/ha_cen), np.log10(oiii_cen/hb_cen), 'ro', alpha = 0.5, markersize = 5)#, label = 'Definite Star Forming')
ax2.set_xlabel(r"$\rm \log([SII]/H\alpha)$", fontsize = 22)
#ax2.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)

#OI/OIII plot
#fig ,ax3= plt.subplots(323)#'OI Scatter Plot')
ax3.plot(refoiha[refoiha < -0.7], o1hamain(refoiha[refoiha < -0.7]),
                  'k', label = 'Ke01 Theoretical Maximum Starburst Line')
ax3.plot(refoiha[refoiha < -0.7], o1hamain(refoiha[refoiha < -0.7]),
                  'k-.', label = 'Ka03 Composite Line')
#ax3.set_ylim(ylims)
ax3.plot(refoiha[refoiha > -1.13], o1halinseyf(refoiha[refoiha > -1.13]),'k--', 
         label = 'Ke06 Liner/Seyfert Division Line')
ax3.set_xlabel(r"$\rm \log([OI]/H\alpha)$", fontsize = 22)
#ax3.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
#ax3.plot(o1ha, o3hb, 'k.', alpha = 0.5, 
#                    markersize = 5, label = 'SF')
#ax3.plot(np.log10(oi_cen/ha_cen), np.log10(oiii_cen/hb_cen), 'ro', alpha = 0.5, markersize = 5)#, label = 'Definite Star Forming')
if error:
    ax1.errorbar(n2ha.flatten(), o3hb.flatten(), xerr = n2ha_err.flatten(),
                yerr = o3hb_err.flatten(),fmt = 'None', marker = 'None', 
                alpha = 0.5, mew = 0, label = 'SF-to-AGN',
                ecolor = 'k', zorder=0)
    ax2.errorbar(s2ha.flatten(), o3hb.flatten(), xerr = s2ha_err.flatten(),
                yerr = o3hb_err.flatten(), fmt = 'None', marker = 'None', c = r,
                alpha = 0.5, mew = 0, label = 'SF-to-AGN', ecolor = 'k',
                zorder=0)
    ax3.errorbar(o1ha.flatten(), o3hb.flatten(), xerr = o1ha_err.flatten(),
                yerr = o3hb_err.flatten(), fmt = 'None', marker = 'None', c = r,
                alpha = 0.5, mew = 0, label = 'SF-to-AGN', ecolor = 'k', 
                zorder=0)
cax = ax1.scatter(n2ha, o3hb,c = r, cmap = cmap)
cax = ax2.scatter(s2ha, o3hb,c = r, cmap = cmap)
cax = ax3.scatter(o1ha, o3hb,c = r, cmap = cmap)
fig.colorbar(cax,extend = 'min')
ax1.plot(np.log10(nii_cen/ha_cen), np.log10(oiii_cen/hb_cen), 'ko', \
         mfc = 'none', mew = 2, markersize = 15)#, label = 'Definite Star Forming')
ax2.plot(np.log10(sii_cen/ha_cen), np.log10(oiii_cen/hb_cen), 'ko', \
         mfc = 'none', mew = 2, markersize = 15)#, label = 'Definite Star Forming')
ax3.plot(np.log10(oi_cen/ha_cen), np.log10(oiii_cen/hb_cen), 'ko', \
         mfc = 'none', mew = 2, markersize = 15)#, label = 'Definite Star Forming')
ax1.scatter(np.log10(nii_cen/ha_cen), np.log10(oiii_cen/hb_cen),marker = 'o', 
            s = 175, c = [0.0], cmap = cmap)#, label = 'Definite Star Forming')
ax2.scatter(np.log10(sii_cen/ha_cen), np.log10(oiii_cen/hb_cen), marker = 'o', \
         s = 175, c = [0.0], cmap = cmap)#, label = 'Definite Star Forming')
ax3.scatter(np.log10(oi_cen/ha_cen), np.log10(oiii_cen/hb_cen), marker = 'o', \
         s = 175, c = [0.0], cmap = cmap)#, label = 'Definite Star Forming')

#ax3.errorbar(-0.9, 0.3, xerr = ratioerror(oi_cen,oi_cen/10,ha_cen,ha_cen/20),
#            yerr = ratioerror(oiii_cen,oiii_cen/20,hb_cen,hb_cen/10), 
#                marker = 'o', c = 'k')
#rs1105
df = pd.read_csv('C:\Users\mugdhapolimera\github\SDSS_spectra\RESOLVE_full_snr5.csv')
df.index = df.name
sdss_oi = df.oi_6300_flux.loc[galname]
sdss_oi_err = df.oi_6300_flux_err.loc[galname]
sdss_oiii =  df.oiii_5007_flux.loc[galname]
sdss_oiii_err = df.oiii_5007_flux_err.loc[galname]
sdss_ha = df.h_alpha_flux.loc[galname]
sdss_ha_err = df.h_alpha_flux_err.loc[galname]
sdss_hb = df.h_beta_flux.loc[galname]
sdss_hb_err = df.h_beta_flux_err.loc[galname]

ax3.errorbar(-1.05, 0.3, xerr = ratioerror(oi_cen,oi_cen_err,ha_cen,ha_cen_err),
                yerr = ratioerror(oiii_cen,oiii_cen_err,hb_cen,hb_cen_err), 
                marker = 'o', c = 'k')
ax3.errorbar(-1.2, 0.3, xerr = ratioerror(sdss_oi,sdss_oi_err,sdss_ha,sdss_ha_err),
            yerr = ratioerror(sdss_oiii,sdss_oiii_err,sdss_hb,sdss_hb_err), 
                marker = 'o', c = 'k')
ax3.errorbar(-0.9, 0.3, xerr = ratioerror(sdss_oi,sdss_oi/23,sdss_ha,sdss_ha/43),
            yerr = ratioerror(sdss_oiii,sdss_oiii/30,sdss_hb,sdss_hb/28), 
                marker = 'o', c = 'k')
ax3.set_xlim(-1.6,-0.8)
ax2.set_xlim(-0.5,-0.1)
ax1.set_xlim(xlims)
ax1.set_ylim(ylims)

#df = pd.read_csv('C:\Users\mugdhapolimera\github\SDSS_spectra\RESOLVE_full_snr5_port.csv')
#df.index = df.name
#galname = 'rs1375'
#sdss_oi = df.Flux_OI_6300.loc[galname]
#sdss_oi_err = df.Flux_OI_6300_Err.loc[galname]
#sdss_sii = df.Flux_SII_6716.loc[galname]+df.Flux_SII_6730.loc[galname]
#sdss_sii_err = df.Flux_SII_6730_Err.loc[galname]
#sdss_oiii =  df.Flux_OIII_5006.loc[galname]
#sdss_oiii_err = df.Flux_OIII_5006_Err.loc[galname]
#sdss_ha = df.Flux_Ha_6562.loc[galname]
#sdss_ha_err = df.Flux_Ha_6562_Err.loc[galname]
#sdss_hb = df.Flux_Hb_4861.loc[galname]
#sdss_hb_err = df.Flux_Hb_4861_Err.loc[galname]
#
#ax3.errorbar(np.log10(sdss_oi/sdss_ha), np.log10(sdss_oiii/sdss_hb), 
#             xerr = ratioerror(sdss_oi,sdss_oi_err,sdss_ha,sdss_ha_err),
#            yerr = ratioerror(sdss_oiii,sdss_oiii_err,sdss_hb,sdss_hb_err), 
#                marker = 'o', c = 'k')
#ax2.errorbar(np.log10(sdss_sii/sdss_ha), np.log10(sdss_oiii/sdss_hb), 
#             xerr = ratioerror(sdss_sii,sdss_sii_err,sdss_ha,sdss_ha_err),
#            yerr = ratioerror(sdss_oiii,sdss_oiii_err,sdss_hb,sdss_hb_err), 
#                marker = 'o', c = 'k')

#ax3.errorbar(np.log10(sdss_oi/sdss_ha),np.log10(sdss_oiii/sdss_hb),
#             xerr = ratioerror(sdss_oi,sdss_oi_err,sdss_ha,sdss_ha_err),
#            yerr = ratioerror(sdss_oiii,sdss_oiii_err,sdss_hb,sdss_hb_err), 
#                marker = 'o', c = 'k')

#df = pd.read_csv(folder+'_smcdext.csv')
##create line ratios/H-alpha and [OIII]/H-beta
#nii_2 = df['nii_6584_flux']
#nii_sum = df['nii_6584_flux']
#nii_sum_err = df['nii_6584_flux_err']**2
## note that the ratio uses only the stronger line, but for S/N reasons we add
## the weaker and multiply by 3/4 since Chris Richardson says the canonical
## line ratio is 3:1 (this needs to be updated with a more precise number)
#oiii = df['oiii_5007_flux']
#oiii_err = df['oiii_5007_flux_err']
#h_alpha = df['h_alpha_flux']
#h_alpha_err = df['h_alpha_flux_err']
#h_beta = df['h_beta_flux']
#h_beta_err = df['h_beta_flux_err']
#oi = df['oi_6300_flux']
#oi_err = df['oi_6300_flux_err']
#sii_sum = df['sii_6717_flux'] + df['sii_6731_flux']
#sii_sum_err = np.sqrt(df['sii_6717_flux_err']**2 + df['sii_6731_flux_err']**2)
#
##Filter Data: all non-negative SEL fluxes and errors; Hbeta >3sigma
#gooddata = ((h_alpha > 0) & (nii_sum > 0) & (oiii > 0) & (oi > 0) &
#            (sii_sum > 0) & (h_beta > 0) & (h_beta > 3*h_beta_err) &
#            (h_alpha_err > 0) & (nii_sum_err > 0) & (oiii_err > 0) & 
#            (oi_err > 0) & (sii_sum_err > 0))
#
#data = gooddata #use ALL galaxy data within catalog
#nii = nii[data]
#nii_sum = nii_sum[data]
#oiii = oiii[data]
#oiii_err = oiii_err[data]
#oi = oi[data]
#oi_err = oi_err[data]
#sii_sum = sii_sum[data]
#sii_sum_err = sii_sum_err[data]
#h_beta = h_beta[data]
#h_beta_err = h_beta_err[data]
#h_alpha = h_alpha[data]
#h_alpha_err = h_alpha_err[data]
#
#fig = plt.figure('NII Scatter Plot - deextincted')
#ax1 = fig.add_subplot(111)
#ax1.set_xlim(-1.5,0.5)
#ax1.set_ylim(-1.0,1.0)
#ax1.plot(refn2ha, n2hamain(refn2ha), 'k', 
#                  label = 'ke01 Theoretical Maximum Starburst Line')
#ax1.plot(refn2ha[refn2ha < 0], n2hacompmin(refn2ha[refn2ha < 0]),
#                      'k-.', label = 'Ka03 Composite Line')
#ax1.plot(n2ha, o3hb, 'k.', alpha = 0.5, markersize = 5)#, label = 'Definite Star Forming')
#ax1.set_xlabel(r"$\rm \log([NII]/H\alpha)$", fontsize = 22)
#ax1.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
#if error: 
#    ax1.errorbar(n2ha.flatten(), o3hb.flatten(), xerr = n2ha_err.flatten(),
#                            yerr = o3hb_err.flatten(), fmt = 'b.', alpha = 0.5,
#                        markersize = 8, mew = 0, label = 'SF-to-AGN', ecolor = 'k')
#
##SII/OIII plot
#fig = plt.figure('SII Scatter Plot - deextincted')
#ax2 = fig.add_subplot(111)
#ax2.plot(refsiiha, s2hamain(refsiiha), 'k',  label = 'Ke01 Line')
#ax2.plot(refsiiha[refsiiha > -0.31], s2halinseyf(refsiiha[refsiiha > -0.31]),
#                  'k--', label = 'Liner/Seyfert Division')
#ax2.set_xlim(-1.5, 0.5)
#ax2.set_ylim(-1.0,1.0)
#ax2.set_xlabel(r"$\rm \log([SII]/H\alpha)$", fontsize = 22)
#ax2.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
#ax2.plot(s2ha, o3hb, 'k.', markersize = 5, 
#                    alpha = 0.5, label = 'SF')
#if error:
#    ax2.errorbar(s2ha.flatten(), o3hb.flatten(), xerr = s2ha_err.flatten(),
#                            yerr = o3hb_err.flatten(), fmt = 'b.', alpha = 0.5,
#                        markersize = 8, mew = 0, label = 'SF-to-AGN', ecolor = 'k')
#
##OI/OIII plot
#fig = plt.figure('OI Scatter Plot - deextincted')
#ax3 = fig.add_subplot(111)
#ax3.plot(refoiha[refoiha < -0.7], o1hamain(refoiha[refoiha < -0.7]),
#                  'k', label = 'Ke01 Theoretical Maximum Starburst Line')
#ax3.plot(refoiha[refoiha < -0.7], o1hamain(refoiha[refoiha < -0.7]),
#                  'k-.', label = 'Ka03 Composite Line')
##ax3.set_xlim(-2.0, -0.4)
##ax3.set_ylim(-1.0,1.0)
#ax3.plot(refoiha[refoiha > -1.13], o1halinseyf(refoiha[refoiha > -1.13]),
#                               'k--', label = 'Ke06 Liner/Seyfert Division Line')
#ax3.set_xlabel(r"$\rm \log([OI]/H\alpha)$", fontsize = 22)
#ax3.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
#ax3.plot(o1ha, o3hb, 'k.', alpha = 0.5, 
#                    markersize = 5, label = 'SF')
#if error:
#    ax3.errorbar(o1ha.flatten(), o3hb.flatten(), xerr = o1ha_err.flatten(),
#                                yerr = o3hb_err.flatten(), fmt = 'b.', alpha = 0.5,
#                            markersize = 8, mew = 0, label = 'SF-to-AGN', ecolor = 'k')
