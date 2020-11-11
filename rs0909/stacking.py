# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 16:33:10 2020

@author: mugdhapolimera
"""

import numpy as np
import pandas as pd
from astropy.table import Table as table
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import patches
import os
from scipy.io import readsav
import itertools
from astropy import wcs
from astropy.coordinates import SkyCoord as sky
import matplotlib as mpl

image1 = fits.open("slin_rs0909_20200523.fits")
image2 = fits.open("slin_rs0909_20200625.fits")
image3 = fits.open("slin_rs0909_20200716.fits")

im1 = image1[0].data
im2 = image2[0].data
im3 = image3[0].data

hdu1 = image1[0].header
hdu2 = image2[0].header
hdu3 = image3[0].header

im1_cen = 213
im2_cen = 204
im3_cen = 267

plt.figure()
plt.plot(im1[213])#hdu1['GALROW']])
plt.plot(im2[204])#hdu2['GALROW']])
plt.plot(im3[267])#hdu3['GALROW']])
plt.plot(im1[213]+im2[204]+im3[267])#im1[hdu1['GALROW']]+im2[hdu2['GALROW']]+im3[hdu3['GALROW']])

im1_new = im1[im1_cen-100:im1_cen+101]
im2_new = im2[im2_cen-100:im2_cen+101]
im3_new = im3[im3_cen-100:im3_cen+101]

plt.figure()
plt.plot(im1_new[100])#hdu1['GALROW']])
plt.plot(im2_new[100])#hdu2['GALROW']])
plt.plot(im3_new[100])#hdu3['GALROW']])
plt.plot(im1_new[100]+im2_new[100]+im3_new[100])#im1[hdu1['GALROW']]+im2[hdu2['GALROW']]+im3[hdu3['GALROW']])

newhdu = fits.PrimaryHDU(data=im1_new,header=hdu1)
newhdu.writeto('slin_rs0909_20200523_crop.fits')

newhdu = fits.PrimaryHDU(data=im2_new,header=hdu2)
newhdu.writeto('slin_rs0909_20200625_crop.fits')

newhdu = fits.PrimaryHDU(data=im3_new,header=hdu3)
newhdu.writeto('slin_rs0909_20200716_crop.fits')