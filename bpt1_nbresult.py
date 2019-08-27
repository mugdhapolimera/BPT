# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 23:09:31 2018

@author: mugdhapolimera
"""

import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt

path = os.path.dirname(os.path.realpath(__file__))
he2_flag = 0
if he2_flag:
    flags = pd.read_csv('resolve_emlineclass_filtered_he2.csv')
else:
    flags = pd.read_csv('C:/Users/mugdhapolimera/github/SDSS_spectra/resolve_emlineclass_filter_new.csv')
inputfile = 'C:/Users/mugdhapolimera/github/SDSS_spectra/RESOLVE_filter_new.csv'
full_df = pd.read_csv(inputfile)
full_df.index = full_df.name
df = full_df.loc[flags.galname]

os.chdir('C:/Users/mugdhapolimera/github/nebulabayes/')
results_nb = pd.read_csv("res_bpass_nicholls_csf/RESOLVE_param_estimates.csv")
Z_index = (results_nb['Parameter'] == 'LOGZ')
full_Z_nb = results_nb[Z_index]
full_Z_nb.index = full_Z_nb['Galaxy Name']
Z_nb = full_Z_nb.loc[flags.galname]

keys = ['defagn', 'composite', 'defstarform']
if he2_flag:
    keys.append('heiisel')
    
#flags['defagn'] = flags['defseyf'] | flags['defliner'] | flags['ambigagn']
marker = {'agntosf': 'g^', 'ambigagn': 'rs', 'composite': 'bs', 'defagn': 'rs', 
          'defliner': 'yo', 'defseyf': 'co', 'heiisel': 'ks',
          'defstarform': 'ko', 'sftoagn': 'm^'}
plt.figure()
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
            
        Z_sel = Z_nb.loc[sel.NAME]['Estimate']
        if len(Z_sel) > 1:
            plt.plot(sel.logmstar, Z_sel+8.76, marker[key], label = key, alpha = 0.5)
            plt.legend(borderaxespad=0., loc = 4)
            plt.xlabel('Stellar Mass')
            plt.ylabel('Metallicity')
            plt.title('M-Z Relation (NebulaBayes + Levesque Grid + All Lines w/o H-alpha)')
m_z = np.loadtxt('C:\Users\mugdhapolimera\github\BPT\M-Z_Tremonti04.txt')
plt.plot(m_z[:,0], m_z[:,1],'r')
#From Manucci 2010 - Polynomial of M-Z relationship marginalized over SFR
m = np.linspace(8.5 - 10, max(df.logmstar) - 10, 100)
z = 8.96 + 0.31*m - 0.23*(m**2) - 0.017*(m**3) + 0.046*(m**4)
#plt.ylim(-0.2, 1.2)
plt.plot(m+10,z)
plt.ylim(7.5,9.1)
