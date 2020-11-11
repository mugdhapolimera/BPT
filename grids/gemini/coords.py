# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 13:07:16 2019

@author: mugdhapolimera

Create Coordinate List text file
"""

import pandas as pd
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import os

os.chdir('C:/Users/mugdhapolimera/github/SDSS_Spectra/')
inputfile = 'RESOLVE_full_blend_dext_new.csv'
resdata = pd.read_csv(inputfile)

res = pd.DataFrame({})
c = SkyCoord(ra=resdata.radeg*u.degree, dec=resdata.dedeg*u.degree)
res['h'] = c.ra.hms.h
res['m'] = c.ra.hms.m
res['s'] = c.ra.hms.s

res['d'] = c.dec.dms.d
res['dm'] = c.dec.dms.m
res['ds'] = c.dec.dms.s

res['search'] = 7*np.ones(len(res))
res['width'] = np.zeros(len(res))

res['name'] = resdata['name']
res.to_csv('RESOLVEfullcoords.txt', sep = " ", index = False)