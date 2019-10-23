# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 14:37:05 2019

@author: mugdhapolimera

Creating Target list for GEMINI FT Proposal October 2019
"""

import numpy as np
import pandas as pd
import os
from scipy.io import readsav

os.chdir('C:\Users\mugdhapolimera\github\SDSS_spectra')

jhu = pd.read_csv('jhu_sfagn.csv', index_col = 'name')
port = pd.read_csv('port_sfagn.csv', index_col = 'name')
nsa = pd.read_csv('nsa_sfagn.csv', index_col = 'name')

unique = np.unique(list(jhu.index) + list(nsa.index) + list(port.index))

print('List of Unique SFing-AGN from JHU, Portsmouth and NSA Catalogs')
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
print(unq)
ra_beg = 22.5
ra_end = 17.5
observable = ((unq.h > ra_beg) & (unq.h < 24)) | (unq.h > 0) & (unq.h < ra_end)

print('\n\nList of SFing-AGN observable between December 1st and Februray 29th')
print('RA between '+ str(ra_beg) + ' and ' + str(ra_end))
print(unq[observable])

target = unq[observable]
resolve = readsav('resolvecatalog.dat')
mu_r = pd.DataFrame(dict(zip(resolve.name, resolve.ifusb)))