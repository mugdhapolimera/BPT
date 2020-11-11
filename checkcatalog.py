# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 07:57:30 2019

@author: mugdhapolimera
"""

import numpy as np
import pandas as pd
from scipy.io import readsav
import os
#os.chdir('C:\Users\mugdhapolimera\github\SDSS_spectra  ')
os.chdir('../SDSS_spectra/')
sfagn = np.array(['rf0117', 'rf0192', 'rf0238', 'rf0338', 'rf0342', 'rf0376',
       'rf0432', 'rf0477', 'rf0503', 'rs0010', 'rs0063', 'rs0070',
       'rs0105', 'rs0111', 'rs0124', 'rs0150', 'rs0421', 'rs0472',
       'rs0545', 'rs0626', 'rs0775', 'rs0909', 'rs1038', 'rs1047',
       'rs1073', 'rs1105', 'rs1143', 'rs1195', 'rs1283', 'rs1292'])


df = readsav('resolvecatalog.dat')

ndx = [x for x in range(len(df.name)) if df.name[x] .decode("utf-8") in targetnames]
for i in ndx:
    print(df.name[i].decode("utf-8"),df.broaddate[i].decode("utf-8"))
    
print('\n Galaxies with no Broad setup data: ')
print(df.name[ndx][np.where(df.broaddate[ndx] == b'NA')[0]])