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
from matplotlib import patches
import matplotlib as mpl
mpl.rcParams.update({'font.size': 20})
mpl.rcParams.update({'axes.linewidth': 2})
mpl.rcParams.update({'lines.linewidth': 2})
#mpl.rcParams['xtick.labelsize'] = label_size 
#mpl.rcParams['ytick.labelsize'] = label_size 
#mpl.rcParams['xtick.length'] = 8
#mpl.rcParams['ytick.length'] = 8
import os
from scipy.io import readsav
import itertools
from astropy import wcs
from astropy.coordinates import SkyCoord as sky
label_size = 15

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
resphot = readsav('../SDSS_spectra/resolvecatalogphot.dat')
df = pd.read_csv('C:\Users\mugdhapolimera\github\SDSS_spectra\RESOLVE_snr5_master.csv')
df.index = df.name
galname = 'rs0010'
#galname = 'rs0022'
#galname = 'rs0756'
if galname == 'rs0775':
    folder = '71146' #rs0775
if galname == 'rs0010':
    folder = '372320' #rs0010
if galname == 'rs0013': 
    folder = '372374' #not SNR 5
if galname == 'rs0756':
    folder = '85507'
if galname == 'rs0022':
    folder = '210660'

os.chdir(r'F:\mugdhapolimera\Documents\UNC\Courses\Research\SAMI Data\/'+ folder)
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
sii_1 = fits.open(folder+'_SII6716_'+filename)[0].data#.flatten()
sii_1_err = fits.open(folder+'_SII6716_'+filename)[1].data#.flatten()
sii_2 = fits.open(folder+'_SII6731_'+filename)[0].data#.flatten()
sii_2_err = fits.open(folder+'_SII6731_'+filename)[1].data#.flatten()
sii = sii_1 + sii_2
velfield = fits.open(folder+'_gas-velocity_default_recom-comp.fits')[0].data#.flatten()
velhdu = fits.open(folder+'_gas-velocity_default_recom-comp.fits')[0].header

cube_hdu = fits.open(folder+'_cube_red.fits')[0].header
cube = fits.open(folder+'_cube_red.fits')[0].data
#cube_hdu = fits.open(folder+'_adaptive_red.fits')[0].header
#cube = fits.open(folder+'_adaptive_red.fits')[0].data

bluecube_hdu = fits.open(folder+'_cube_blue.fits')[0].header
bluecube = fits.open(folder+'_cube_blue.fits')[0].data#.flatten()
#bluecube_hdu = fits.open(folder+'_adaptive_blue.fits')[0].header
#bluecube = fits.open(folder+'_adaptive_blue.fits')[0].data#.flatten()
pixelscale = cube_hdu['CDELT1'] #deg/pix
v = resdata['vhel'][np.where(resdata['name'] == galname)] #km/s
z = v/3e5
d = (v/70)*10**6 #Hubble's law - d in pc
pc = 2*np.pi*d*(pixelscale/360)/1000 #kpc/pix
b_a = np.float(df.loc[galname].b_a)
#pc = pixelscale *3600 #in arcsec/pix
MAX = np.max(halpha[~np.isnan(halpha)])
center = [[int(cube_hdu['CRPIX1'])],[int(cube_hdu['CRPIX2'])]]#np.where(halpha == MAX)
#center = [center[0][0], center[1][0]]
dx = np.arange(np.shape(halpha)[0])- center[0]
dy = np.arange(np.shape(halpha)[1])- center[1]
dxy = np.array(zip(dx,dy))
r = np.zeros(halpha.shape)
for i in range(len(dx)):
    for j in range(len(dy)):
        r[i][j] = np.sqrt((pc*dx[i])**2 + (pc*-dy[j]/b_a)**2)

hdu = cube_hdu
lam0 = cube_hdu['CRVAL3']-((cube_hdu['CRPIX3']-1)*cube_hdu['CDELT3'])
lam = (np.arange(cube_hdu['NAXIS3']))*cube_hdu['CDELT3'] + lam0
bluelam0 = bluecube_hdu['CRVAL3']-((bluecube_hdu['CRPIX3']-1)*bluecube_hdu['CDELT3'])
bluelam = (np.arange(bluecube_hdu['NAXIS3']))*bluecube_hdu['CDELT3'] + bluelam0
###############################################################################
# Calculate S/N from the continuum of the blue spectra
###############################################################################
start = 4500
end = 4800
startndx = np.where(abs(bluelam-start) == min(abs(bluelam-start)))[0][0]
endndx = np.where(abs(bluelam-end) == min(abs(bluelam-end)))[0][0]
bluesnr = np.zeros(np.shape(bluecube[0,:,:]))
blue = fits.open(folder+'_cube_blue.fits')
bluevar = blue[1].data
for i,j in itertools.product(list(range(50)),list(range(50))):
    spaxel = [i,j]
    signal = np.nanmedian(bluecube[startndx:endndx,spaxel[0],spaxel[1]]) 
    #in units of 10**(-16) erg /s /cm**2 /angstrom /pixel
    noise = np.sqrt(np.nanmedian(bluevar[startndx:endndx,spaxel[0],spaxel[1]]))
    bluesnr[i][j] = signal/noise

start = 6475
end = 6615
startndx = np.where(abs(lam-start) == min(abs(lam-start)))[0][0]
endndx = np.where(abs(lam-end) == min(abs(lam-end)))[0][0]
redsnr = np.zeros(np.shape(halpha))
spaxel = [25,25] #center
red = fits.open(folder+'_cube_red.fits')
redvar = red[1].data
for i,j in itertools.product(list(range(50)),list(range(50))):
    spaxel = [i,j]
    signal = 2*np.nanmedian(cube[startndx:endndx,spaxel[0],spaxel[1]]) 
    #in units of 10**(-16) erg /s /cm**2 /angstrom /pixel
    noise = np.sqrt(2*np.nanmedian(bluevar[startndx:endndx,spaxel[0],spaxel[1]]))
    redsnr[i][j] = signal/noise

#plt.figure()
#plt.plot(lam,cube[:,25,25])
#plt.plot(bluelam,bluecube[:,25,25])

#Try to create "fake" SDSS spectrum from SAMI spectrum
#rndx = np.where(r<=2.05)
#sdss = cube[67:1990,rndx[0],rndx[1]]
#nans = np.where(np.isnan(sdss))
#sdss[nans] = 0
#sdss2arc = np.sum(sdss,axis=1)
#sdss2arc = sdss2arc.T[:,None,None]
#hdu['NAXIS'] = 3
#hdu['NAXIS1'] = sdss2arc.shape[1]
#hdu['NAXIS2'] = sdss2arc.shape[2]
#hdu['NAXIS3'] = sdss2arc.shape[0]
#
##hdu['CRVAL3'] = hdu['CRVAL1']
#hdu['CRVAL1'] = 1
#hdu['CRVAL2'] = 1
#
#hdu['CRPIX3'] = 1
#hdu['CRVAL3'] = lam[67]
#
##hdu['CDELT3'] = hdu['CDELT1']
##hdu['CD3_3'] = hdu['CD1_1']
#hdu['CDELT1'] = 1
#hdu['CDELT2'] = 1


mesh = np.meshgrid(np.arange(50), np.arange(50))
#To search for a particular wavelength
line = (1+z)*6561 #6300
ndx = np.where(abs(lam-line) == min(abs(lam-line)))[0][0]

ha_cen = halpha[center]
hb_cen = hbeta[center]
nii_cen = nii[center]
sii_cen = sii[center]
sii_cen_err = np.sqrt(sii_1[center]**2 + sii_2[center]**2)
oi_cen = oi[center]
oiii_cen = oiii[center]
ha_cen_err = halpha_err[center]
hb_cen_err = hbeta_err[center]
nii_cen_err = nii_err[center]
oi_cen_err = oi_err[center]
oiii_cen_err = oiii_err[center]
#center = [center[0][0],center[1][0]]

nans = [~np.isnan(halpha) & ~np.isnan(hbeta) & ~np.isnan(nii) & ~np.isnan(sii_1)
        & ~np.isnan(sii_2) & ~np.isnan(oi) & ~np.isnan(oiii)]
snr = 5 #np.nanmedian(oi/oi_err)
#err = [(nii/nii_err > snr) & (sii_1/sii_1_err > snr) & (sii_2/sii_2_err > snr) & 
#       (oi/oi_err > snr) & (oiii/oiii_err > snr) & (hbeta/hbeta_err > snr) & 
#       (halpha/halpha_err > snr)]
err = [(halpha/halpha_err > snr)]
snr = 10
contsnr = (bluesnr > snr) & (redsnr > snr)
good = nans and err and contsnr
#good = good[0]
ndxs = np.array(list(itertools.product(np.arange(50),np.arange(50))))
ndxs = ndxs[good.ravel()]
#ax.scatter(ndxs[:,0], ndxs[:,1], color = 'white', marker = 's')
badndxs = np.array(list(itertools.product(np.arange(50),np.arange(50))))
badndx = badndxs[~good.ravel()]

from copy import copy
cmap = copy(plt.cm.seismic)
cmap.set_bad('gray',0.8)

w = wcs.WCS(cube_hdu)
w = w.dropaxis(2)
xy = np.meshgrid(np.arange(50), np.arange(50))
coords = w.all_pix2world(xy[0], xy[1],0)

image = np.ma.masked_where(~good,np.log10(oi/halpha))
#image = np.ma.masked_where(~good,(halpha))

goodimage = np.ma.masked_where(((bluesnr < 5) & (redsnr < 5)),np.log10(oi/halpha))
#goodimage = np.ma.masked_where(((bluesnr < 5) & (redsnr < 5)),(halpha))
goodimage[:,0:15] = np.nan
goodimage[:,34:] = np.nan
fig = plt.figure()
ax = plt.subplot(projection = w)
#ax.imshow(cube[ndx,:,:],norm = colors.Normalize(vmin = 0, vmax = 0.04439), 
#          cmap = 'inferno')
cax = ax.imshow(goodimage,
                norm = colors.Normalize(vmin = -1.4, vmax = -1.0), 
                cmap = 'seismic')
cax = ax.imshow(image,norm = colors.Normalize(vmin = -1.4, vmax = -1.0), 
          cmap = cmap)
#cax = ax.imshow(goodimage,
#                norm = colors.Normalize(vmin = 0.09, vmax = 1.0), 
#                cmap = 'seismic')
#cax = ax.imshow(image,norm = colors.Normalize(vmin = 0.09, vmax = 1.0), 
#          cmap = cmap)
ax.plot(25,25,'x', c = 'black')
#cmap = colors.ListedColormap(['purple','orange','blue','yellow','green','red',\
#                              'cyan','black','white','lime','purple'])#,\
#                              'magenta','magenta','magenta','magenta','magenta',
#                              'magenta','magenta','magenta','magenta','magenta'])
#cmap = plt.get_cmap('rainbow_r',18)#int(np.max(r)))
#boundaries = np.array([0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])#,
#                       5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0])
#norm = colors.BoundaryNorm(boundaries, cmap.N, clip=True)
#cax = ax.scatter(mesh[0],mesh[1], s = 3, c = r, cmap = cmap, norm = norm)

#Want to plot annuli at every 0.5 kpc; 5 pixels = 1 kpc
#for i in range(0,30,5):
#    #To add annulus with 0.5 kpc diameter 
##    e1 = patches.Ellipse((25,25), width = b_a*0.5*i, height = 0.5*i, 
##                         linewidth = 2, fill = False, zorder = 2)
##    ax.add_patch(e1)
#    #To add annulus with 1 kpc radius
#    e1 = patches.Ellipse((25,25), width = b_a*i*2, height = 2*i, 
#                         linewidth = 2, fill = False, zorder = 2, color = 'g')
#    ax.add_patch(e1)
Re = resphot.radr50p[resdata.name == galname] #in arcsec
Re_deg = Re/3600.0
height = 2*Re/(pixelscale*3600)
Re_kpc = 2*np.pi*d*(Re_deg/360)/1000 #kpc/pix

e_Re = patches.Ellipse((25,25), width = b_a*height, height = height, 
                         linewidth = 2, fill = False, zorder = 2, ls = '--',
                         color = 'black')
ax.add_patch(e_Re)
e1 = patches.Ellipse((25,25), width = b_a*5*2, height = 2*5, 
                     linewidth = 2, fill = False, zorder = 2, color = 'yellow')
ax.add_patch(e1)
    
cb = fig.colorbar(cax,extend = 'min')
ax.set_xlabel('RA')
ax.set_ylabel('Dec')
ax.coords[0].set_ticklabel(size=10.0)
ax.coords[1].set_ticklabel(size=10.0)
cb.ax.tick_params(labelsize=10.0)
#fig.xticks(fontsize=15)
#ax.tick_params(axis = 'both',labelsize = 10)
#ax.scatter(badndx[:,1], badndx[:,0], color = 'gray', marker = 's', alpha = 0.7,
#           edgecolors = 'none')

goodvel = np.ma.masked_where(((bluesnr < 5) & (redsnr < 5)),(velfield[0,:,:]))

fig = plt.figure()
ax = plt.subplot(projection = w)
cax = ax.imshow(velfield[0,:,:],
                norm = colors.Normalize(vmin = -75, vmax = 75), 
                cmap = 'seismic')
fig.colorbar(cax,extend = 'min')
ax.set_xlabel('RA')
ax.set_ylabel('Dec')
ax.plot(25,25,'x', c = 'black')

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
fig,(ax1,ax2,ax3) = plt.subplots(1,3,sharey = True)

#plt.title(galname+' SNR > '+str(snr))
good = (r <= 3.05)
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

#catid = [folder]*len(halpha)
#data = list(zip(catid, hbeta, hbeta_err,oiii, oiii_err, oi, oi_err,
#                halpha, halpha_err, nii, nii_err, sii_1, sii_1_err, sii_2, sii_2_err))
#names = ['CATID',
#         'h_beta_flux', 'h_beta_flux_err', 
#       'oiii_5007_flux', 'oiii_5007_flux_err',
#       'oi_6300_flux', 'oi_6300_flux_err', 
#       'h_alpha_flux','h_alpha_flux_err',
#       'nii_6584_flux', 'nii_6584_flux_err', 
#       'sii_6717_flux','sii_6717_flux_err',
#       'sii_6731_flux', 'sii_6731_flux_err']
#
#df = pd.DataFrame(data, columns = names)
#df.to_csv(folder+'.csv')
sii = sii_1 + sii_2
sii_err = np.sqrt(sii_1_err**2 + sii_2_err**2)
n2ha = np.log10(nii/halpha)
s2ha = np.log10(sii/halpha)
o1ha= np.log10(oi/halpha)
o3hb = np.log10(oiii/hbeta)
n2ha_err = ratioerror(nii, nii_err, halpha, halpha_err)
s2ha_err = ratioerror(sii, sii_err, halpha, halpha_err)
o1ha_err = ratioerror(oi, oi_err, halpha, halpha_err)
o3hb_err = ratioerror(oiii, oiii_err, hbeta, hbeta_err)

refn2ha = np.linspace(-3.0, 0.35)
refoiha = np.linspace(-2.5, -0.4)
refsiiha = np.linspace(-2, 0.32,100)

xlims = [-1.0,0.0]
ylims = [-0.15,0.4]
cmap = plt.get_cmap('rainbow_r',10)#int(np.max(r)/0.5))
boundaries = np.array([0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5])#, 4.0])#, 4.5, 5.0])
norm = colors.BoundaryNorm(boundaries, cmap.N, clip=True)

#NII plot
ax1.plot(refn2ha, n2hamain(refn2ha), 'k', 
                  label = 'ke01 Theoretical Maximum Starburst Line')
ax1.plot(refn2ha[refn2ha < 0], n2hacompmin(refn2ha[refn2ha < 0]),
                      'k-.', label = 'Ka03 Composite Line')
ax1.set_xlabel(r"$\rm \log([NII]/H\alpha)$", fontsize = 22)
ax1.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)

#SII plot
ax2.plot(refsiiha, s2hamain(refsiiha), 'k',  label = 'Ke01 Line')
ax2.plot(refsiiha[refsiiha > -0.32], s2halinseyf(refsiiha[refsiiha > -0.32]),
                  'k--', label = 'Liner/Seyfert Division')
ax2.set_xlabel(r"$\rm \log([SII]/H\alpha)$", fontsize = 22)
#ax2.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)

#OI plot
ax3.plot(refoiha[refoiha < -0.7], o1hamain(refoiha[refoiha < -0.7]),
                  'k', label = 'Ke01 Theoretical Maximum Starburst Line')
ax3.plot(refoiha[refoiha < -0.7], o1hamain(refoiha[refoiha < -0.7]),
                  'k-.', label = 'Ka03 Composite Line')
ax3.plot(refoiha[refoiha > -1.13], o1halinseyf(refoiha[refoiha > -1.13]),'k--', 
         label = 'Ke06 Liner/Seyfert Division Line')
ax3.set_xlabel(r"$\rm \log([OI]/H\alpha)$", fontsize = 22)
#ax3.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)

if error:
    ax1.errorbar(n2ha.flatten(), o3hb.flatten(), xerr = n2ha_err.flatten(),
                yerr = o3hb_err.flatten(),fmt = 'None', marker = 'None', 
                alpha = 0.5, mew = 0, label = 'SF-to-AGN',
                ecolor = 'k', zorder=0, norm = norm)
    ax2.errorbar(s2ha.flatten(), o3hb.flatten(), xerr = s2ha_err.flatten(),
                yerr = o3hb_err.flatten(), fmt = 'None', marker = 'None', c = r,
                alpha = 0.5, mew = 0, label = 'SF-to-AGN', ecolor = 'k',
                zorder=0, norm = norm)
    ax3.errorbar(o1ha.flatten(), o3hb.flatten(), xerr = o1ha_err.flatten(),
                yerr = o3hb_err.flatten(), fmt = 'None', marker = 'None', c = r,
                alpha = 0.5, mew = 0, label = 'SF-to-AGN', ecolor = 'k', 
                zorder=0, norm = norm)
cax = ax1.scatter(n2ha, o3hb,c = r, cmap = cmap, norm = norm)
cax = ax2.scatter(s2ha, o3hb,c = r, cmap = cmap, norm = norm)
cax = ax3.scatter(o1ha, o3hb,c = r, cmap = cmap, norm = norm)
fig.colorbar(cax,extend = 'min', ticks = boundaries, boundaries= boundaries)

#Plotting the center with a bigger marker
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

df = pd.read_csv('C:\Users\mugdhapolimera\github\SDSS_spectra\RESOLVE_snr5_master.csv')
df.index = df.name
sdss_oi = df.oi_6300_flux.loc[galname]
sdss_oi_err = df.oi_6300_flux_err.loc[galname]
sdss_oiii =  df.oiii_5007_flux.loc[galname]
sdss_oiii_err = df.oiii_5007_flux_err.loc[galname]
sdss_ha = df.h_alpha_flux.loc[galname]
sdss_ha_err = df.h_alpha_flux_err.loc[galname]
sdss_hb = df.h_beta_flux.loc[galname]
sdss_hb_err = df.h_beta_flux_err.loc[galname]
ax1.tick_params(length=8, width=2)
ax2.tick_params(length=8, width=2)
ax3.tick_params(length=8, width=2)

#Plotting SAMI errorbar
#ax3.errorbar(-1.05, 0.3, xerr = ratioerror(oi_cen,oi_cen_err,ha_cen,ha_cen_err),
#                yerr = ratioerror(oiii_cen,oiii_cen_err,hb_cen,hb_cen_err), 
#                marker = 'o', c = 'k')
##Plotting SDSS errorbar
#ax3.errorbar(-1.2, 0.3, xerr = ratioerror(sdss_oi,sdss_oi_err,sdss_ha,sdss_ha_err),
#            yerr = ratioerror(sdss_oiii,sdss_oiii_err,sdss_hb,sdss_hb_err), 
#                marker = 'o', c = 'k')
##Plotting expected GEMINI errorbar
#ax3.errorbar(-0.9, 0.3, xerr = ratioerror(sdss_oi,sdss_oi/37,sdss_ha,sdss_ha/78),
#            yerr = ratioerror(sdss_oiii,sdss_oiii/40,sdss_hb,sdss_hb/33), 
#                marker = 'o', c = 'k')
ax3.set_xlim(-1.6,-0.8)
ax2.set_xlim(-0.5,-0.1)
ax1.set_xlim(xlims)
ax1.set_ylim(ylims)

#line = (1+z)*6561 #6300
#ndx = np.where(abs(lam-line) == min(abs(lam-line)))[0][0]
#spec = fits.open('372320_annular_red.fits')[0].data
#goodspec = np.ma.masked_where(np.isnan(spec[ndx,:,:]),spec[ndx,:,:])
#
#cmap = copy(plt.cm.seismic)
#cmap.set_bad('gray',0.8)
#
#fig = plt.figure()
#ax = plt.subplot(projection = w)
#cax = ax.imshow(goodspec,
#                norm = colors.Normalize(vmin = np.nanmin(goodspec), vmax = np.nanmax(goodspec)), 
#                cmap = cmap)
#fig.colorbar(cax,extend = 'min')
#ax.set_xlabel('RA')
#ax.set_ylabel('Dec')
#ax.plot(25,25,'x', c = 'black')
#plt.figure()
#plt.plot(lam,spec[:,25,25])
#
#spec = fits.open('372320_spectrum_1-4-arcsec_red.fits')[0].data
#plt.figure()
#plt.plot(lam/(1+z),spec*(1+z),'k')
#plt.axvline(x = 6562.8,ymin = 0, ymax = 1,color = 'orange',ls = '--')
#plt.axvline(x = 6583.45,ymin = 0, ymax = 1,color = 'orange',ls = '--')
#plt.axvline(x = 6548.05,ymin = 0, ymax = 1,color = 'orange',ls = '--')
#plt.xlabel('Wavelength (Angstroms)')
#plt.ylabel('Flux (arbitraty units)')
#
#spec = fits.open('372320_spectrum_re_blue.fits')[0].data
#plt.figure()
#plt.plot(bluelam,spec)