# -*- coding: utf-8 -*-
"""
Created on Mon Dec 03 14:02:38 2018

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
    flags = pd.read_csv('resolve_emlineclass_bpt1.csv')
inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_SDSS_filtered.pkl'
full_df = pd.read_pickle(inputfile)
df = full_df.loc[flags.galname]

os.chdir('C:/Users/mugdhapolimera/github/nebulabayes/')
results_nb = pd.read_csv("results_izi_prior/RESOLVE_param_estimates.csv")
Z_index = (results_nb['Parameter'] == 'LOGZ')
full_Z_nb = results_nb[Z_index]
full_Z_nb.index = full_Z_nb['Galaxy Name']
Z_nb = full_Z_nb.loc[flags.galname]

keys = ['agntosf', 'defagn', 'composite', 'defstarform', 'sftoagn']
if he2_flag:
    keys.append('heiisel')
    
flags['defagn'] = flags['defseyf'] | flags['defliner'] | flags['ambigagn']
marker = {'agntosf': 'g^', 'ambigagn': 'rs', 'composite': 'bs', 'defagn': 'rs', 
          'defliner': 'yo', 'defseyf': 'co', 'heiisel': 'ks',
          'defstarform': 'k.', 'sftoagn': 'm^'}
plt.figure()
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
            
        Z_sel = Z_nb.loc[sel.NAME]['Estimate']
        
        plt.plot(sel.logmstar, Z_sel, marker[key], label = key)
        plt.legend(borderaxespad=0., loc = 4)
        plt.xlabel('Stellar Mass')
        plt.ylabel('Metallicity')
        plt.title('M-Z Relation using NebulaBayes + Levesque Grid')

m_z = np.loadtxt('C:\Users\mugdhapolimera\github\BPT\BPT\BPT\M-Z_Tremonti04.txt')
plt.plot(m_z[:,0], m_z[:,1],'r')
#From Manucci 2010 - Polynomial of M-Z relationship marginalized over SFR
m = np.linspace(8.5 - 10, max(df.logmstar) - 10, 100)
z = 8.96 + 0.31*m - 0.23*(m**2) - 0.017*(m**3) + 0.046*(m**4)
#plt.ylim(-0.2, 1.2)
plt.plot(m+10,z)

os.chdir('C:/Anaconda2/Lib/site-packages/NebulaBayes/docs/')
results_nb = pd.read_csv("results_bpass_full/RESOLVE_param_estimates.csv")
Z_index = (results_nb['Parameter'] == 'LOGZ')
full_Z_nb = results_nb[Z_index]
full_Z_nb.index = full_Z_nb['Galaxy Name']
Z_nb = full_Z_nb.loc[flags.galname]

keys = ['agntosf', 'defagn', 'composite', 'defstarform', 'sftoagn']
if he2_flag:
    keys.append('heiisel')
    
flags['defagn'] = flags['defseyf'] | flags['defliner'] | flags['ambigagn']
marker = {'agntosf': 'g^', 'ambigagn': 'rs', 'composite': 'bs', 'defagn': 'rs', 
          'defliner': 'yo', 'defseyf': 'co', 'heiisel': 'ks',
          'defstarform': 'k.', 'sftoagn': 'm^'}
plt.figure()
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
            
        Z_sel = Z_nb.loc[sel.NAME]['Estimate']
        
        plt.plot(sel.logmstar, Z_sel, marker[key], label = key)
        plt.legend(borderaxespad=0., loc = 4)
        plt.xlabel('Stellar Mass')
        plt.ylabel('Metallicity')
        plt.title('M-Z Relation using NebulaBayes + BPASS')

m_z = np.loadtxt('C:\Users\mugdhapolimera\github\BPT\BPT\BPT\M-Z_Tremonti04.txt')
plt.plot(m_z[:,0], m_z[:,1],'r')
#From Manucci 2010 - Polynomial of M-Z relationship marginalized over SFR
m = np.linspace(8.5 - 10, max(df.logmstar) - 10, 100)
z = 8.96 + 0.31*m - 0.23*(m**2) - 0.017*(m**3) + 0.046*(m**4)
#plt.ylim(-0.2, 1.2)
plt.plot(m+10,z)

os.chdir('C:/Anaconda2/Lib/site-packages/NebulaBayes/docs/')
results_nb = pd.read_csv("results_k13/RESOLVE_param_estimates.csv")
Z_index = (results_nb['Parameter'] == 'LOGZ')
full_Z_nb = results_nb[Z_index]
full_Z_nb.index = full_Z_nb['Galaxy Name']
Z_nb = full_Z_nb.loc[flags.galname]

keys = ['agntosf', 'defagn', 'composite', 'defstarform', 'sftoagn']
if he2_flag:
    keys.append('heiisel')
    
flags['defagn'] = flags['defseyf'] | flags['defliner'] | flags['ambigagn']
marker = {'agntosf': 'g^', 'ambigagn': 'rs', 'composite': 'bs', 'defagn': 'rs', 
          'defliner': 'yo', 'defseyf': 'co', 'heiisel': 'ks',
          'defstarform': 'k.', 'sftoagn': 'm^'}
plt.figure()
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
            
        Z_sel = Z_nb.loc[sel.NAME]['Estimate']
        
        plt.plot(sel.logmstar, Z_sel, marker[key], label = key)
        plt.legend(borderaxespad=0., loc = 4)
        plt.xlabel('Stellar Mass')
        plt.ylabel('Metallicity')
        plt.title('M-Z Relation using NebulaBayes + K13 Grid')

m_z = np.loadtxt('C:\Users\mugdhapolimera\github\BPT\BPT\BPT\M-Z_Tremonti04.txt')
plt.plot(m_z[:,0], m_z[:,1],'r')
#From Manucci 2010 - Polynomial of M-Z relationship marginalized over SFR
m = np.linspace(8.5 - 10, max(df.logmstar) - 10, 100)
z = 8.96 + 0.31*m - 0.23*(m**2) - 0.017*(m**3) + 0.046*(m**4)
#plt.ylim(-0.2, 1.2)
plt.plot(m+10,z)

os.chdir('C:/Anaconda2/Lib/site-packages/NebulaBayes/docs/')
results_nb = pd.read_csv("results_k13_deredden/RESOLVE_param_estimates.csv")
Z_index = (results_nb['Parameter'] == 'LOGZ')
full_Z_nb = results_nb[Z_index]
full_Z_nb.index = full_Z_nb['Galaxy Name']
Z_nb = full_Z_nb.loc[flags.galname]

keys = ['agntosf', 'defagn', 'composite', 'defstarform', 'sftoagn']
if he2_flag:
    keys.append('heiisel')
    
flags['defagn'] = flags['defseyf'] | flags['defliner'] | flags['ambigagn']
marker = {'agntosf': 'g^', 'ambigagn': 'rs', 'composite': 'bs', 'defagn': 'rs', 
          'defliner': 'yo', 'defseyf': 'co', 'heiisel': 'ks',
          'defstarform': 'k.', 'sftoagn': 'm^'}
plt.figure()
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
            
        Z_sel = Z_nb.loc[sel.NAME]['Estimate']
        
        plt.plot(sel.logmstar, Z_sel, marker[key], label = key)
        plt.legend(borderaxespad=0., loc = 4)
        plt.xlabel('Stellar Mass')
        plt.ylabel('Metallicity')
        plt.title('M-Z Relation using NebulaBayes + K13 Grid (NB Deredden with Prior)')

m_z = np.loadtxt('C:\Users\mugdhapolimera\github\BPT\BPT\BPT\M-Z_Tremonti04.txt')
plt.plot(m_z[:,0], m_z[:,1],'r')
#From Manucci 2010 - Polynomial of M-Z relationship marginalized over SFR
m = np.linspace(8.5 - 10, max(df.logmstar) - 10, 100)
z = 8.96 + 0.31*m - 0.23*(m**2) - 0.017*(m**3) + 0.046*(m**4)
#plt.ylim(-0.2, 1.2)
plt.plot(m+10,z)

os.chdir('C:/Anaconda2/Lib/site-packages/NebulaBayes/docs/')
results_nb = pd.read_csv("results_k13_prior/RESOLVE_param_estimates.csv")
Z_index = (results_nb['Parameter'] == 'LOGZ')
full_Z_nb = results_nb[Z_index]
full_Z_nb.index = full_Z_nb['Galaxy Name']
Z_nb = full_Z_nb.loc[flags.galname]

keys = ['agntosf', 'defagn', 'composite', 'defstarform', 'sftoagn']
if he2_flag:
    keys.append('heiisel')
    
flags['defagn'] = flags['defseyf'] | flags['defliner'] | flags['ambigagn']
marker = {'agntosf': 'g^', 'ambigagn': 'rs', 'composite': 'bs', 'defagn': 'rs', 
          'defliner': 'yo', 'defseyf': 'co', 'heiisel': 'ks',
          'defstarform': 'k.', 'sftoagn': 'm^'}
plt.figure()
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
            
        Z_sel = Z_nb.loc[sel.NAME]['Estimate']
        
        plt.plot(sel.logmstar, Z_sel, marker[key], label = key)
        plt.legend(borderaxespad=0., loc = 4)
        plt.xlabel('Stellar Mass')
        plt.ylabel('Metallicity')
        plt.title('M-Z Relation using NebulaBayes + K13 Grid (with Prior)')

m_z = np.loadtxt('C:\Users\mugdhapolimera\github\BPT\BPT\BPT\M-Z_Tremonti04.txt')
plt.plot(m_z[:,0], m_z[:,1],'r')
#From Manucci 2010 - Polynomial of M-Z relationship marginalized over SFR
m = np.linspace(8.5 - 10, max(df.logmstar) - 10, 100)
z = 8.96 + 0.31*m - 0.23*(m**2) - 0.017*(m**3) + 0.046*(m**4)
#plt.ylim(-0.2, 1.2)
plt.plot(m+10,z)

os.chdir('C:/Anaconda2/Lib/site-packages/NebulaBayes/docs/')
results_nb = pd.read_csv("results_nbgrid/RESOLVE_param_estimates.csv")
Z_index = (results_nb['Parameter'] == 'LOGZ')
full_Z_nb = results_nb[Z_index]
full_Z_nb.index = full_Z_nb['Galaxy Name']
Z_nb = full_Z_nb.loc[flags.galname]

keys = ['agntosf', 'defagn', 'composite', 'defstarform', 'sftoagn']
if he2_flag:
    keys.append('heiisel')
    
flags['defagn'] = flags['defseyf'] | flags['defliner'] | flags['ambigagn']
marker = {'agntosf': 'g^', 'ambigagn': 'rs', 'composite': 'bs', 'defagn': 'rs', 
          'defliner': 'yo', 'defseyf': 'co', 'heiisel': 'ks',
          'defstarform': 'k.', 'sftoagn': 'm^'}
plt.figure()
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
            
        Z_sel = Z_nb.loc[sel.NAME]['Estimate']
        
        plt.plot(sel.logmstar, Z_sel, marker[key], label = key)
        plt.legend(borderaxespad=0., loc = 4)
        plt.xlabel('Stellar Mass')
        plt.ylabel('Metallicity')
        plt.title('M-Z Relation using NebulaBayes + NB MAPPINGS V Grid')

m_z = np.loadtxt('C:\Users\mugdhapolimera\github\BPT\BPT\BPT\M-Z_Tremonti04.txt')
plt.plot(m_z[:,0], m_z[:,1],'r')
#From Manucci 2010 - Polynomial of M-Z relationship marginalized over SFR
m = np.linspace(8.5 - 10, max(df.logmstar) - 10, 100)
z = 8.96 + 0.31*m - 0.23*(m**2) - 0.017*(m**3) + 0.046*(m**4)
#plt.ylim(-0.2, 1.2)
plt.plot(m+10,z)

os.chdir('C:/Anaconda2/Lib/site-packages/NebulaBayes/docs/')
results_nb = pd.read_csv("results_nbgrid_deredden/RESOLVE_param_estimates.csv")
Z_index = (results_nb['Parameter'] == 'LOGZ')
full_Z_nb = results_nb[Z_index]
full_Z_nb.index = full_Z_nb['Galaxy Name']
Z_nb = full_Z_nb.loc[flags.galname]

keys = ['agntosf', 'defagn', 'composite', 'defstarform', 'sftoagn']
if he2_flag:
    keys.append('heiisel')
    
flags['defagn'] = flags['defseyf'] | flags['defliner'] | flags['ambigagn']
marker = {'agntosf': 'g^', 'ambigagn': 'rs', 'composite': 'bs', 'defagn': 'rs', 
          'defliner': 'yo', 'defseyf': 'co', 'heiisel': 'ks',
          'defstarform': 'k.', 'sftoagn': 'm^'}
plt.figure()
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
            
        Z_sel = Z_nb.loc[sel.NAME]['Estimate']
        
        plt.plot(sel.logmstar, Z_sel, marker[key], label = key)
        plt.legend(borderaxespad=0., loc = 4)
        plt.xlabel('Stellar Mass')
        plt.ylabel('Metallicity')
        plt.title('M-Z Relation using NebulaBayes + NB MAPPINGS V Grid (NB Deredden with Prior)')

m_z = np.loadtxt('C:\Users\mugdhapolimera\github\BPT\BPT\BPT\M-Z_Tremonti04.txt')
plt.plot(m_z[:,0], m_z[:,1],'r')
#From Manucci 2010 - Polynomial of M-Z relationship marginalized over SFR
m = np.linspace(8.5 - 10, max(df.logmstar) - 10, 100)
z = 8.96 + 0.31*m - 0.23*(m**2) - 0.017*(m**3) + 0.046*(m**4)
#plt.ylim(-0.2, 1.2)
plt.plot(m+10,z)

os.chdir('C:/Anaconda2/Lib/site-packages/NebulaBayes/docs/')
results_nb = pd.read_csv("results_nbgrid_prior/RESOLVE_param_estimates.csv")
Z_index = (results_nb['Parameter'] == 'LOGZ')
full_Z_nb = results_nb[Z_index]
full_Z_nb.index = full_Z_nb['Galaxy Name']
Z_nb = full_Z_nb.loc[flags.galname]

keys = ['agntosf', 'defagn', 'composite', 'defstarform', 'sftoagn']
if he2_flag:
    keys.append('heiisel')
    
flags['defagn'] = flags['defseyf'] | flags['defliner'] | flags['ambigagn']
marker = {'agntosf': 'g^', 'ambigagn': 'rs', 'composite': 'bs', 'defagn': 'rs', 
          'defliner': 'yo', 'defseyf': 'co', 'heiisel': 'ks',
          'defstarform': 'k.', 'sftoagn': 'm^'}
plt.figure()
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
            
        Z_sel = Z_nb.loc[sel.NAME]['Estimate']
        
        plt.plot(sel.logmstar, Z_sel, marker[key], label = key)
        plt.legend(borderaxespad=0., loc = 4)
        plt.xlabel('Stellar Mass')
        plt.ylabel('Metallicity')
        plt.title('M-Z Relation using NebulaBayes + NB MAPPINGS V Grid (with Prior)')

m_z = np.loadtxt('C:\Users\mugdhapolimera\github\BPT\BPT\BPT\M-Z_Tremonti04.txt')
plt.plot(m_z[:,0], m_z[:,1],'r')
#From Manucci 2010 - Polynomial of M-Z relationship marginalized over SFR
m = np.linspace(8.5 - 10, max(df.logmstar) - 10, 100)
z = 8.96 + 0.31*m - 0.23*(m**2) - 0.017*(m**3) + 0.046*(m**4)
#plt.ylim(-0.2, 1.2)
plt.plot(m+10,z)

os.chdir('C:/Anaconda2/Lib/site-packages/NebulaBayes/docs/')
results_nb = pd.read_csv("results_izi_deredden/RESOLVE_param_estimates.csv")
Z_index = (results_nb['Parameter'] == 'LOGZ')
full_Z_nb = results_nb[Z_index]
full_Z_nb.index = full_Z_nb['Galaxy Name']
Z_nb = full_Z_nb.loc[flags.galname]

keys = ['agntosf', 'defagn', 'composite', 'defstarform', 'sftoagn']
if he2_flag:
    keys.append('heiisel')
    
flags['defagn'] = flags['defseyf'] | flags['defliner'] | flags['ambigagn']
marker = {'agntosf': 'g^', 'ambigagn': 'rs', 'composite': 'bs', 'defagn': 'rs', 
          'defliner': 'yo', 'defseyf': 'co', 'heiisel': 'ks',
          'defstarform': 'k.', 'sftoagn': 'm^'}
plt.figure()
for key in keys:
        sel = df.iloc[np.where(flags[key])[0]]
            
        Z_sel = Z_nb.loc[sel.NAME]['Estimate']
        
        plt.plot(sel.logmstar, Z_sel, marker[key], label = key)
        plt.legend(borderaxespad=0., loc = 4)
        plt.xlabel('Stellar Mass')
        plt.ylabel('Metallicity')
        plt.title('M-Z Relation using NebulaBayes + Levesque + Deredden')

m_z = np.loadtxt('C:\Users\mugdhapolimera\github\BPT\BPT\BPT\M-Z_Tremonti04.txt')
plt.plot(m_z[:,0], m_z[:,1],'r')
#From Manucci 2010 - Polynomial of M-Z relationship marginalized over SFR
m = np.linspace(8.5 - 10, max(df.logmstar) - 10, 100)
z = 8.96 + 0.31*m - 0.23*(m**2) - 0.017*(m**3) + 0.046*(m**4)
#plt.ylim(-0.2, 1.2)
plt.plot(m+10,z)
