
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 11:13:54 2019

@author: mugdhapolimera

Creating target list for GEMINI FT proposal

"""
import pandas as pd
import os
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy.io.idl import readsav
import numpy as np

os.chdir('C:\Users\mugdhapolimera\github\SDSS_spectra/')
resolve = pd.read_csv('RESOLVE_full_raw.csv', index_col = 'name')

jhuflag = pd.read_csv('resolve_emlineclass_full_snr5_jhu.csv', index_col = 'galname')
portflag = pd.read_csv('resolve_emlineclass_full_snr5_port.csv', index_col = 'galname')
nsaflag = pd.read_csv('resolve_emlineclass_full_snr5_nsa.csv', index_col = 'galname')

port = pd.read_csv('RESOLVE_full_snr5_port.csv', index_col = 'name')[portflag.sftoagn]
jhu = pd.read_csv('RESOLVE_full_snr5.csv', index_col = 'name')[jhuflag.sftoagn]
nsa = pd.read_csv('NSA_RESOLVE.csv', index_col = 'resname').loc[nsaflag.index.values[nsaflag.sftoagn]]
 
port = port[['radeg','dedeg']]   
jhu = jhu[['radeg','dedeg']]   
nsa = nsa[['radeg','dedeg']]

allunq = np.unique(list(jhuflag.index) + list(nsaflag.index) + list(portflag.index))   
#sftoagn = df[ambigsel1 & dwarf][['radeg','dedeg']]

#sfagn = pd.read_csv('uniquesfagn.csv')
unique = np.unique(list(jhu.index) + list(nsa.index) + list(port.index))

#print('List of Unique SFing-AGN from JHU, Portsmouth and NSA Catalogs')
#print(len(unique))
#print (jhu.loc[[x for x in jhu.index if x in unique]])
#print (port.loc[[x for x in port.index \
#                 if (x in unique) & (x not in jhu.index)]])
#print (nsa.loc[[x for x in nsa.index \
#                if (x in unique) & (x not in jhu.index) & (x not in port.index)]])

unq = jhu.loc[[x for x in jhu.index if x in unique]]
unq = unq.append(port.loc[[x for x in port.index \
                 if (x in unique) & (x not in jhu.index)]])
unq = unq.append(nsa.loc[[x for x in nsa.index \
                if (x in unique) & (x not in jhu.index) & (x not in port.index)]])

sfagn = unq

c = SkyCoord(ra=sfagn.radeg*u.degree, dec=sfagn.dedeg*u.degree)
sfagn['h'] = c.ra.hms.h
sfagn['m'] = c.ra.hms.m
sfagn['s'] = c.ra.hms.s

sfagn['d'] = c.dec.dms.d
sfagn['dm'] = c.dec.dms.m
sfagn['ds'] = c.dec.dms.s

#print(sfagn)

#target range 
#November 2019 deadline- obs from jan to march 
ra_low = 8.9
ra_high = 18.5
sel = (ra_low < sfagn.h) & (sfagn.h< ra_high)

#print(sfagn[sel])
targetnames = list(sfagn[sel].index.values)
#resolve = readsav('../../../SDSS_spectra/resolvecatalog.dat')
#mu_r = pd.DataFrame(dict(zip(resolve.name, resolve.ifusb)))

ra = []
dec = []
rsys = []
for x in nov_target:
    ra.append(str(int(sfagn.loc[x].h))+':'+\
              str(int(sfagn.loc[x].m))+':'+\
              str(sfagn.loc[x].s))
    dec.append(str(int(sfagn.loc[x].d))+':'+\
              str(int(sfagn.loc[x].dm))+':'+\
              str(sfagn.loc[x].ds))
    rsys.append('Vega')
    #ndx = np.where(resdata.name == x)
    #mu_r = resdata.ifusb[ndx]
M_r = list(df.loc[nov_target]['absrmag'])                
data = zip(nov_target,ra,dec, M_r,rsys)
nov = pd.DataFrame(data = data, columns = ['Name','RAJ200','DecJ200','r','r_sys'])
nov.to_csv('November_targetlist.csv')