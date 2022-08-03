# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 07:19:45 2021

@author: mugdhapolimera
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord as sky
from matplotlib import colors
from scipy.io import readsav
from matplotlib import patches
from astropy import units as u
from astropy.visualization.wcsaxes import SphericalCircle
from astropy.coordinates import Angle
import os
resdb = readsav('../SDSS_spectra/resolvecatalog.dat')
resphotdb = readsav('../SDSS_spectra/resolvecatalogphot.dat')
ecodb = readsav('../SDSS_spectra/eco_wresa_032918.dat')

galname = 'rs0966'
#galname = 'rs1139'
origfile = galname+'-w1-int_old.fits'
newfile = galname+'-w1-int_new.fits'

origfile = 'mask_'+galname+'-w1-int_old.fits'
newfile = 'mask_'+galname+'-w1-int_new.fits'
##
#origfile = galname+'-wise-w1-modelpsfh_wisemask_old.fit'
#newfile = galname+'-wise-w1-modelpsfh_wisemask_new.fit'


os.chdir('../SDSS_spectra/mid_ir/')
resra = resdb.ra[resdb.name == galname]
resdec = resdb.dec[resdb.name == galname]
resradius = Angle(resphotdb.radr24p5[resdb.name == galname]*2, unit = 'arcsec')
#resaprad = Angle(130.0688* 0.40697825554, unit = 'arcsec') #old apread with unnecessary factor of 2
#resaprad = Angle(65.034425* 0.40697825554, unit = 'arcsec')
resaprad = resradius
res = fits.open(origfile)[0].data
reshdu = fits.open(origfile)[0].header
#res = fits.open('swarptest_bgsub_crop_rs0075_1421p015_ac51-w1-int-3.fits')[0].data
#reshdu = fits.open('swarptest_bgsub_crop_rs0075_1421p015_ac51-w1-int-3.fits')[0].header
#res = fits.open('magzptestbgsub_crop_rs0075_1421p015_ac51-w1-int-3.fits')[0].data
#reshdu = fits.open('magzptestbgsub_crop_rs0075_1421p015_ac51-w1-int-3.fits')[0].header
#res = fits.open('mask_rs0075-w1-int.fits')[0].data
#reshdu = fits.open('mask_rs0075-w1-int.fits')[0].header
#res = -2.5*np.log10(res)+reshdu['MAGZP']

wres = wcs.WCS(reshdu)
fig = plt.figure()
ax1 = plt.subplot(projection = wres)
cax = ax1.imshow(res, cmap = 'seismic')#,
#                 norm = colors.Normalize(vmin = 0, vmax = 12)) 
ax1.plot(resra,resdec, 'x', color = 'white', transform=ax1.get_transform('world'))
e1 = SphericalCircle((resra,resdec)*u.deg, resradius, 
                         linewidth = 2, fill = False, zorder = 2, color = 'g',
                         transform=ax1.get_transform('fk5'))
ax1.add_patch(e1)
e2 = SphericalCircle((resra,resdec)*u.deg,resaprad, 
                         linewidth = 2, fill = False, zorder = 2, color = 'yellow',
                         transform=ax1.get_transform('fk5'))
ax1.add_patch(e2)
#ax1.set_xlim(110,200)
#ax1.set_ylim(110,200)

eco= fits.open(newfile)[0].data
ecohdu = fits.open(newfile)[0].header
              
#eco = fits.open('bgsub_crop_rs0075_1421p015_ac51-w1-int-3.fits')[0].data
#ecohdu = fits.open('bgsub_crop_rs0075_1421p015_ac51-w1-int-3.fits')[0].header
#eco = fits.open('magzptestmask_rs0075-w1-int.fits')[0].data
#ecohdu = fits.open('magzptestmask_rs0075-w1-int.fits')[0].header
#eco= -2.5*np.log10(eco)+ecohdu['MAGZP']

weco = wcs.WCS(ecohdu)
fig = plt.figure()
ax1 = plt.subplot(projection = weco)
cax = ax1.imshow(eco, cmap = 'seismic')#, 
#                 norm = colors.Normalize(vmin = 0, vmax = 12)) 
ax1.plot(resra,resdec, 'x', color = 'white', transform=ax1.get_transform('world'))
#ax1.set_xlim(125,175)
#ax1.set_ylim(125,175)
e1 = SphericalCircle((resra,resdec)*u.deg, resradius, 
                         linewidth = 2, fill = False, zorder = 2, color = 'g',
                         transform=ax1.get_transform('fk5'))
ax1.add_patch(e1)
e2 = SphericalCircle((resra,resdec)*u.deg,resaprad, 
                         linewidth = 2, fill = False, zorder = 2, color = 'yellow',
                         transform=ax1.get_transform('fk5'))
ax1.add_patch(e2)
#ax1.set_xlim(110,200)
#ax1.set_ylim(110,200)

fig = plt.figure()
ax1 = plt.subplot(projection = weco)
cax = ax1.imshow(eco-res, cmap = 'seismic', 
                 norm = colors.Normalize(vmin = 0, vmax = 1)) 
ax1.plot(resra,resdec, 'x', color = 'white', transform=ax1.get_transform('world'))
#ax1.set_xlim(125,175)
#ax1.set_ylim(125,175)
e1 = SphericalCircle((resra,resdec)*u.deg, resradius, 
                         linewidth = 2, fill = False, zorder = 2, color = 'g',
                         transform=ax1.get_transform('fk5'))
ax1.add_patch(e1)
e2 = SphericalCircle((resra,resdec)*u.deg,resaprad, 
                         linewidth = 2, fill = False, zorder = 2, color = 'yellow',
                         transform=ax1.get_transform('fk5'))
ax1.add_patch(e2)
#ax1.set_xlim(110,200)
#ax1.set_ylim(110,200)
