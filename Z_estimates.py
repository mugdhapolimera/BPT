# -*- coding: utf-8 -*-
"""
Created on Mon Oct 01 13:44:53 2018

@author: mugdhapolimera
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
pd.set_option('display.max_columns', 500)
bpt = pd.read_csv("C:/Users/mugdhapolimera/github/BPT/BPT/BPT/resolve_emlineclass_new.csv")

#IZI + NebulaBayes
izi = pd.read_csv("C:/Anaconda2/Lib/site-packages/NebulaBayes/docs/results_izi/RESOLVE_param_estimates_izi.csv")
Z_index = (izi['Parameter'] == 'LOGZ')
izi_Z = izi[Z_index]
izi_Z.index = range(len(izi_Z))
Z_izi = izi_Z['Estimate']
izi_Z['CI68_low'] = izi_Z['CI68_low'].replace('#NAME?', '-inf').astype(np.float)
izi_Z['CI68_high'] = izi_Z['CI68_high'].replace('np.inf', 'inf').astype(np.float)

errup_izi = izi_Z['CI68_high'] - izi_Z['Estimate']
errdown_izi = izi_Z['Estimate'] - izi_Z['CI68_low']

print 'HeII Selected'
print izi_Z[bpt['heiisel']]

print 'Ambiguous AGN'
print izi_Z[bpt['ambigagn']]